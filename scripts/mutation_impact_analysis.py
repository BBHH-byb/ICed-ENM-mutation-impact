#!/usr/bin/env python3
"""Run the ICed-ENM mutation impact analysis pipeline."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import math
import os
import shutil
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any


AMINO_TYPES = (
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLU",
    "GLN",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
)


@dataclass(frozen=True)
class Residue:
    chain: str
    resid: str
    resname: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class MutationTask:
    pdb_id: str
    analysis_chains: str
    chain: str
    resid: str
    from_resname: str
    to_resname: str
    work_dir: Path
    mutant_dir: Path
    modeller_path: Path | None
    seed: int
    max_attempts: int


@dataclass(frozen=True)
class ModellerApi:
    Alignment: Any
    Environ: Any
    Model: Any
    Selection: Any
    log: Any
    autosched: Any
    ConjugateGradients: Any
    MolecularDynamics: Any


def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed

def positive_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive real number") from exc
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("must be a positive real number")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run mutation impact analysis using ICed-ENM normal modes, "
            "mutant structure generation, and vibrational entropy."
        )
    )
    parser.add_argument("input_pdb", type=Path, help="Input PDB file, e.g. 1CRN.pdb.")
    parser.add_argument("chain", help="Chain ID(s) to analyze, e.g. A or ABC.")

    entropy_group = parser.add_mutually_exclusive_group()
    entropy_group.add_argument(
        "--entropy-percent",
        type=positive_float,
        default=5.0,
        help="Percentage of CA modes to use for vibrational entropy. Default: 5.",
    )
    entropy_group.add_argument(
        "--entropy-mode",
        type=positive_int,
        help="Absolute number of modes to use for vibrational entropy.",
    )

    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=4,
        help="Maximum number of parallel mutation and NMA jobs. Default: 4.",
    )
    parser.add_argument(
        "--nma-core",
        type=positive_int,
        default=1,
        help=(
            "OpenMP core count passed to each ICed_ENM_NMA run. "
            "Default: 1 to avoid oversubscribing cores during parallel mutant NMA."
        ),
    )
    parser.add_argument(
        "--max-attempts",
        type=positive_int,
        default=3,
        help="Maximum mutation-generation attempts per mutant. Default: 3.",
    )
    parser.add_argument(
        "--seed",
        type=positive_int,
        default=49837,
        help="Base random seed for mutant generation. Default: 49837.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: <pdb_id>_<chain>.",
    )
    parser.add_argument(
        "--iced-enm-bin",
        type=Path,
        default=Path("./ICed_ENM_NMA"),
        help="Path to the ICed_ENM_NMA executable. Default: ./ICed_ENM_NMA.",
    )
    parser.add_argument(
        "--modeller-path",
        type=Path,
        default=None,
        help=(
            "Optional path to the Modeller Python package, site-packages "
            "directory, modlib directory, or Modeller install root. "
            "Default: use Modeller from the active Python environment."
        ),
    )
    return parser.parse_args()

def stable_seed(base_seed: int, *parts: str) -> int:
    key = "|".join((str(base_seed), *parts)).encode()
    digest = hashlib.sha256(key).hexdigest()
    return int(digest[:8], 16) % 900_000 + 1

def run_command(command: list[str], cwd: Path, log_file: Path | None = None) -> None:
    proc = subprocess.run(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if log_file is not None:
        log_file.write_text(proc.stdout)
    if proc.returncode != 0:
        message = f"Command failed with exit code {proc.returncode}: {' '.join(command)}"
        if log_file is not None:
            message += f"\nSee log: {log_file}"
        raise RuntimeError(message)


def eval_dir(output_dir: Path) -> Path:
    path = output_dir / "evals"
    path.mkdir(parents=True, exist_ok=True)
    return path


def logs_dir(output_dir: Path) -> Path:
    path = output_dir / "logs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def eval_file_path(output_dir: Path, prefix: str) -> Path:
    return eval_dir(output_dir) / f"{prefix}_eval.txt"


def eval_meta_path(output_dir: Path, prefix: str) -> Path:
    return eval_dir(output_dir) / f"{prefix}_eval.meta"


def legacy_eval_file_path(output_dir: Path, prefix: str) -> Path:
    return output_dir / f"{prefix}_eval.txt"


def collect_eval_output(output_dir: Path, prefix: str) -> Path:
    destination = eval_file_path(output_dir, prefix)
    legacy = legacy_eval_file_path(output_dir, prefix)
    if legacy.exists():
        if destination.exists():
            destination.unlink()
        shutil.move(str(legacy), destination)
    return destination


def find_or_collect_eval_output(output_dir: Path, prefix: str) -> Path:
    destination = eval_file_path(output_dir, prefix)
    if destination.exists():
        return destination
    return collect_eval_output(output_dir, prefix)


def eval_metadata_matches(output_dir: Path, prefix: str, chains: str, nmodes: int) -> bool:
    meta_file = eval_meta_path(output_dir, prefix)
    if not meta_file.exists():
        return False
    metadata: dict[str, str] = {}
    for line in meta_file.read_text().splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        metadata[key.strip()] = value.strip()
    return metadata.get("chains") == chains and metadata.get("nmodes") == str(nmodes)


def write_eval_metadata(output_dir: Path, prefix: str, chains: str, nmodes: int) -> None:
    eval_meta_path(output_dir, prefix).write_text(f"chains={chains}\nnmodes={nmodes}\n")


def remove_stale_eval(output_dir: Path, prefix: str) -> None:
    eval_file_path(output_dir, prefix).unlink(missing_ok=True)
    legacy_eval_file_path(output_dir, prefix).unlink(missing_ok=True)
    eval_meta_path(output_dir, prefix).unlink(missing_ok=True)


def read_ca_residues(ref_ca_pdb: Path) -> list[Residue]:
    residues: list[Residue] = []
    seen: set[tuple[str, str]] = set()
    with ref_ca_pdb.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16] != " CA ":
                continue
            chain = line[21:22].strip()
            resid = line[22:26].strip()
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            key = (chain, resid)
            if key in seen:
                continue
            seen.add(key)
            residues.append(Residue(chain=chain, resid=resid, resname=resname, x=x, y=y, z=z))
    return residues

def run_wt_iced_enm(args: argparse.Namespace, output_dir: Path, nmodes: int | None = None) -> Path:
    log_file = logs_dir(output_dir) / "wt_ICed_ENM_NMA.log"
    command = [
        str(args.iced_enm_bin.resolve()),
        str(args.input_pdb.resolve()),
        args.chain,
        "--core",
        str(args.nma_core),
    ]
    if nmodes is not None:
        command.extend(["--mode", str(nmodes)])
    command.extend(["--out-prefix", "wt"])
    run_command(command, cwd=output_dir, log_file=log_file)
    collect_eval_output(output_dir, "wt")
    ref_ca = output_dir / "wt_ref_CA.pdb"
    if not ref_ca.exists():
        raise RuntimeError(f"WT ICed-ENM did not generate {ref_ca}")
    return ref_ca


def cleanup_nma_outputs(output_dir: Path, prefix: str, keep_ref_ca: bool = False) -> None:
    removable = [
        output_dir / f"{prefix}_evec_CA.txt",
        output_dir / f"{prefix}_RMSF_CA.txt",
    ]
    if not keep_ref_ca:
        removable.append(output_dir / f"{prefix}_ref_CA.pdb")
    for path in removable:
        path.unlink(missing_ok=True)


def refine_mutation_selection(atmsel, api: ModellerApi) -> None:
    md = api.MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return="FINAL")
    init_vel = True
    schedule = (
        (200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
        (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0)),
    )
    for iterations, equilibration, temperatures in schedule:
        for temperature in temperatures:
            md.optimize(
                atmsel,
                init_velocities=init_vel,
                temperature=temperature,
                max_iterations=iterations,
                equilibrate=equilibration,
            )
            init_vel = False


def optimize_mutation_selection(atmsel, schedule, api: ModellerApi) -> None:
    for step in schedule:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    refine_mutation_selection(atmsel, api)
    cg = api.ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def make_mutation_restraints(model, alignment, api: ModellerApi) -> None:
    restraints = model.restraints
    restraints.clear()
    selection = api.Selection(model)
    for restraint_type in ("stereo", "phi-psi_binormal"):
        restraints.make(selection, restraint_type=restraint_type, aln=alignment, spline_on_site=True)
    for restraint_type in ("omega", "chi1", "chi2", "chi3", "chi4"):
        restraints.make(
            selection,
            restraint_type=f"{restraint_type}_dihedral",
            spline_range=4.0,
            spline_dx=0.3,
            spline_min_points=5,
            aln=alignment,
            spline_on_site=True,
        )


def normalized_modeller_path(modeller_path: Path) -> Path:
    resolved = modeller_path.resolve()
    if (resolved / "modlib" / "modeller").is_dir():
        return resolved / "modlib"
    return resolved


def configure_modeller_path(modeller_path: Path | None) -> None:
    if modeller_path is not None:
        resolved = modeller_path.resolve()
        import_path = normalized_modeller_path(resolved)
        sys.path.insert(0, str(import_path))
        os.environ["PYTHONPATH"] = (
            str(import_path)
            if not os.environ.get("PYTHONPATH")
            else str(import_path) + os.pathsep + os.environ["PYTHONPATH"]
        )
        if (resolved / "modlib" / "modeller").is_dir():
            version_token = resolved.name.split("-")[-1].replace(".", "v")
            os.environ[f"MODINSTALL{version_token}"] = str(resolved) + os.sep
        return

    active_prefix = str(Path(sys.prefix).resolve())
    cleaned_sys_path: list[str] = []
    for entry in sys.path:
        if not entry:
            cleaned_sys_path.append(entry)
            continue
        resolved = str(Path(entry).resolve())
        if "modeller" in resolved.lower() and not resolved.startswith(active_prefix):
            continue
        cleaned_sys_path.append(entry)
    sys.path[:] = cleaned_sys_path

    pythonpath = os.environ.get("PYTHONPATH")
    if pythonpath:
        cleaned_pythonpath = []
        for entry in pythonpath.split(os.pathsep):
            if not entry:
                continue
            resolved = str(Path(entry).resolve())
            if "modeller" in resolved.lower() and not resolved.startswith(active_prefix):
                continue
            cleaned_pythonpath.append(entry)
        if cleaned_pythonpath:
            os.environ["PYTHONPATH"] = os.pathsep.join(cleaned_pythonpath)
        else:
            os.environ.pop("PYTHONPATH", None)


def load_modeller(modeller_path: Path | None) -> ModellerApi:
    configure_modeller_path(modeller_path)
    try:
        from modeller import Alignment, Environ, Model, Selection, log
        from modeller.automodel import autosched
        from modeller.optimizers import ConjugateGradients, MolecularDynamics
    except Exception as exc:
        raise RuntimeError(
            "Failed to import Modeller. Activate the conda environment that has "
            "Modeller installed and make sure the Modeller license key is set."
        ) from exc

    return ModellerApi(
        Alignment=Alignment,
        Environ=Environ,
        Model=Model,
        Selection=Selection,
        log=log,
        autosched=autosched,
        ConjugateGradients=ConjugateGradients,
        MolecularDynamics=MolecularDynamics,
    )


def mutate_model(
    modelname: str,
    respos: str,
    restyp: str,
    chain: str,
    randseed: int,
    modeller_path: Path | None,
) -> None:
    api = load_modeller(modeller_path)

    api.log.verbose()
    env = api.Environ(rand_seed=-randseed)
    env.io.hetatm = True
    env.edat.dynamic_sphere = False
    env.edat.dynamic_lennard = True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")

    model = api.Model(env, file=modelname)
    alignment = api.Alignment(env)
    alignment.append_model(model, atom_files=modelname, align_codes=modelname)

    selection = api.Selection(model.chains[chain].residues[respos])
    selection.mutate(residue_type=restyp)
    alignment.append_model(model, align_codes=modelname)

    model.clear_topology()
    model.generate_topology(alignment[-1])
    model.transfer_xyz(alignment)
    model.build(initialize_xyz=False, build_method="INTERNAL_COORDINATES")

    original_model = api.Model(env, file=modelname)
    model.res_num_from(original_model, alignment)

    tmp_file = f"{modelname}{chain}{restyp}{respos}.tmp"
    output_file = f"{modelname}{chain}{respos}{restyp}.pdb"
    model.write(file=tmp_file)
    model.read(file=tmp_file)

    make_mutation_restraints(model, alignment, api)
    model.env.edat.nonbonded_sel_atoms = 1
    schedule = api.autosched.loop.make_for_model(model)

    selection = api.Selection(model.chains[chain].residues[respos])
    model.restraints.unpick_all()
    model.restraints.pick(selection)

    selection.energy()
    selection.randomize_xyz(deviation=4.0)

    model.env.edat.nonbonded_sel_atoms = 2
    optimize_mutation_selection(selection, schedule, api)

    model.env.edat.nonbonded_sel_atoms = 1
    optimize_mutation_selection(selection, schedule, api)

    selection.energy()
    model.write(file=output_file)
    Path(tmp_file).unlink(missing_ok=True)


def make_mutation_tasks(
    pdb_id: str,
    analysis_chains: str,
    residues: list[Residue],
    output_dir: Path,
    modeller_path: Path | None,
    seed: int,
    max_attempts: int,
) -> list[MutationTask]:
    mutant_dir = output_dir / "mutant_pdbs"
    tasks: list[MutationTask] = []
    for residue in residues:
        for to_resname in AMINO_TYPES:
            if to_resname == residue.resname:
                continue
            task_seed = stable_seed(seed, pdb_id, residue.chain, residue.resid, to_resname)
            tasks.append(
                MutationTask(
                    pdb_id=pdb_id,
                    analysis_chains=analysis_chains,
                    chain=residue.chain,
                    resid=residue.resid,
                    from_resname=residue.resname,
                    to_resname=to_resname,
                    work_dir=output_dir,
                    mutant_dir=mutant_dir,
                    modeller_path=modeller_path,
                    seed=task_seed,
                    max_attempts=max_attempts,
                )
            )
    return tasks

def generate_mutant(task: MutationTask) -> tuple[MutationTask, bool, str]:
    task.mutant_dir.mkdir(parents=True, exist_ok=True)
    target_name = f"{task.pdb_id}{task.chain}{task.resid}{task.to_resname}.pdb"
    target_path = task.mutant_dir / target_name
    local_output = task.work_dir / target_name
    log_file = logs_dir(task.work_dir) / f"{task.pdb_id}{task.chain}{task.resid}{task.to_resname}.log"

    if target_path.exists():
        log_file.write_text(f"Skipped: existing mutant PDB found at {target_path}\n")
        return task, True, str(target_path)

    local_output.unlink(missing_ok=True)

    logs: list[str] = []
    for attempt in range(1, task.max_attempts + 1):
        attempt_seed = task.seed + attempt - 1
        stdout = io.StringIO()
        return_code = 0
        error_text = ""
        current_dir = Path.cwd()
        try:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stdout):
                os.chdir(task.work_dir)
                mutate_model(
                    task.pdb_id,
                    task.resid,
                    task.to_resname,
                    task.chain,
                    attempt_seed,
                    task.modeller_path,
                )
        except Exception:
            return_code = 1
            error_text = traceback.format_exc()
        finally:
            os.chdir(current_dir)
        logs.append(
            f"Attempt {attempt}/{task.max_attempts}: "
            f"{task.chain} {task.resid} {task.from_resname}->{task.to_resname}, "
            f"seed={attempt_seed}, exit={return_code}\n{stdout.getvalue()}{error_text}"
        )
        if local_output.exists():
            shutil.move(str(local_output), target_path)
            log_file.write_text("\n".join(logs) + f"\nSuccess: {target_path}\n")
            return task, True, str(target_path)

    log_file.write_text("\n".join(logs) + "\nFailed.\n")
    return task, False, str(log_file)

def generate_all_mutants(tasks: list[MutationTask], jobs: int, failed_list: Path) -> None:
    failed: list[str] = []
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(generate_mutant, task) for task in tasks]
        for future in as_completed(futures):
            task, ok, detail = future.result()
            if ok:
                print(f"Ready {task.chain} {task.resid} {task.from_resname}->{task.to_resname}: {detail}")
            else:
                message = (
                    f"{task.pdb_id} {task.chain} {task.resid} "
                    f"{task.from_resname} {task.to_resname} {detail}"
                )
                print(f"Failed {message}")
                failed.append(message)
    failed_list.write_text("\n".join(failed) + ("\n" if failed else ""))


def determine_entropy_nmodes(args: argparse.Namespace, n_residues: int) -> int:
    if args.entropy_mode is not None:
        return args.entropy_mode
    nmodes = round(n_residues * args.entropy_percent / 100.0)
    return max(nmodes, 1)


def run_mutant_nma(
    task: MutationTask,
    iced_enm_bin: Path,
    nmodes: int,
    nma_core: int,
) -> tuple[MutationTask, bool, str]:
    prefix = f"{task.chain}{task.resid}{task.to_resname}"
    mutant_pdb = (task.mutant_dir / f"{task.pdb_id}{task.chain}{task.resid}{task.to_resname}.pdb").resolve()
    log_file = logs_dir(task.work_dir) / f"{prefix}_ICed_ENM_NMA.log"

    if not mutant_pdb.exists():
        return task, False, f"missing mutant PDB: {mutant_pdb}"

    expected_eval = find_or_collect_eval_output(task.work_dir, prefix)
    if expected_eval.exists():
        try:
            if len(read_eval_file(expected_eval)) >= nmodes and eval_metadata_matches(
                task.work_dir,
                prefix,
                task.analysis_chains,
                nmodes,
            ):
                return task, True, f"existing NMA output: {expected_eval}"
        except ValueError:
            pass
    remove_stale_eval(task.work_dir, prefix)

    command = [
        str(iced_enm_bin),
        str(mutant_pdb),
        task.analysis_chains,
        "--core",
        str(nma_core),
        "--mode",
        str(nmodes),
        "--out-prefix",
        prefix,
    ]
    try:
        run_command(command, cwd=task.work_dir, log_file=log_file)
    except RuntimeError as exc:
        return task, False, str(exc)

    expected_eval = collect_eval_output(task.work_dir, prefix)
    if not expected_eval.exists():
        return task, False, f"missing NMA output: {expected_eval}"
    write_eval_metadata(task.work_dir, prefix, task.analysis_chains, nmodes)
    cleanup_nma_outputs(task.work_dir, prefix)
    return task, True, str(expected_eval)


def run_all_mutant_nma(
    tasks: list[MutationTask],
    iced_enm_bin: Path,
    nmodes: int,
    nma_core: int,
    jobs: int,
    failed_list: Path,
) -> None:
    failed: list[str] = []
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(run_mutant_nma, task, iced_enm_bin, nmodes, nma_core)
            for task in tasks
        ]
        for future in as_completed(futures):
            task, ok, detail = future.result()
            prefix = f"{task.chain}{task.resid}{task.to_resname}"
            if ok:
                print(f"Calculated NMA {prefix}: {detail}")
            else:
                message = (
                    f"{task.pdb_id} {task.chain} {task.resid} "
                    f"{task.from_resname} {task.to_resname} {detail}"
                )
                print(f"Failed NMA {message}")
                failed.append(message)
    failed_list.write_text("\n".join(failed) + ("\n" if failed else ""))


def read_eval_file(eval_file: Path) -> list[float]:
    values: list[float] = []
    with eval_file.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                values.append(float(stripped.split()[0]))
    return values


def vibrational_entropy_score(wt_eval: list[float], mut_eval: list[float], nmodes: int) -> float:
    if len(wt_eval) < nmodes:
        raise ValueError(f"WT eval has {len(wt_eval)} modes, but {nmodes} are required.")
    if len(mut_eval) < nmodes:
        raise ValueError(f"Mutant eval has {len(mut_eval)} modes, but {nmodes} are required.")

    wt_part = wt_eval[:nmodes]
    mut_part = mut_eval[:nmodes]
    if any(value <= 0.0 for value in wt_part):
        raise ValueError("WT eigenvalues must be positive for vibrational entropy calculation.")
    if any(value <= 0.0 for value in mut_part):
        raise ValueError("Mutant eigenvalues must be positive for vibrational entropy calculation.")

    return sum(math.log(value) for value in wt_part) - sum(math.log(value) for value in mut_part)


def distance(a: Residue, b: Residue) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def residue_valid_mask(residues: list[Residue], contact_cutoff: float = 8.0, degree_min: int = 6) -> list[bool]:
    valid: list[bool] = []
    for i, residue_i in enumerate(residues):
        degree = 0
        for j, residue_j in enumerate(residues):
            if i == j:
                continue
            if distance(residue_i, residue_j) <= contact_cutoff:
                degree += 1
        valid.append(degree >= degree_min)
    return valid


def gaussian_smooth_scores_self_include(
    residues: list[Residue],
    scores: list[float],
    valid_mask: list[bool],
    sigma: float = 5.0,
    radius: float = 15.0,
) -> list[float]:
    if len(residues) != len(scores) or len(scores) != len(valid_mask):
        raise ValueError("residues, scores, and valid_mask lengths must match.")

    smoothed: list[float] = []
    for i, residue_i in enumerate(residues):
        weighted_sum = 0.0
        weight_total = 0.0
        for j, residue_j in enumerate(residues):
            if not valid_mask[j]:
                continue
            if not valid_mask[i] and i == j:
                continue
            dij = distance(residue_i, residue_j)
            if dij > radius:
                continue
            weight = math.exp(-(dij**2) / (2.0 * sigma**2))
            weighted_sum += weight * scores[j]
            weight_total += weight
        smoothed.append(weighted_sum / weight_total if weight_total > 0.0 else math.nan)
    return smoothed


def percentile(values: list[float], q: float) -> float:
    finite_values = sorted(value for value in values if not math.isnan(value))
    if not finite_values:
        raise ValueError("Cannot calculate percentile of an empty score set.")
    if q <= 0.0:
        return finite_values[0]
    if q >= 100.0:
        return finite_values[-1]
    position = (len(finite_values) - 1) * q / 100.0
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return finite_values[int(position)]
    weight = position - lower
    return finite_values[lower] * (1.0 - weight) + finite_values[upper] * weight


def make_visualization_score_map(
    residue_score_map: dict[tuple[str, str], float],
    high_percentile: float = 90.0,
) -> dict[tuple[str, str], float]:
    values = [score for score in residue_score_map.values() if not math.isnan(score)]
    if not values:
        return {}

    low = min(values)
    high = percentile(values, high_percentile)
    if high <= low:
        return {key: 0.0 for key in residue_score_map}

    visualization_scores: dict[tuple[str, str], float] = {}
    for key, score in residue_score_map.items():
        if math.isnan(score):
            visualization_scores[key] = math.nan
            continue
        clipped = min(max(score, low), high)
        visualization_scores[key] = 100.0 * (clipped - low) / (high - low)
    return visualization_scores


def write_residue_impact_bfactor_pdb(
    input_pdb: Path,
    output_pdb: Path,
    residue_score_map: dict[tuple[str, str], float],
    missing_score: float = 0.0,
) -> None:
    with input_pdb.open() as in_handle, output_pdb.open("w") as out_handle:
        for line in in_handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 66:
                chain = line[21:22].strip()
                resid = line[22:26].strip()
                score = residue_score_map.get((chain, resid), missing_score)
                if math.isnan(score):
                    score = missing_score
                out_handle.write(f"{line[:60]}{score:6.2f}{line[66:]}")
            else:
                out_handle.write(line)


def write_mutation_impact_table(
    tasks: list[MutationTask],
    residues: list[Residue],
    output_dir: Path,
    input_pdb: Path,
    wt_eval_file: Path,
    nmodes: int,
    mutation_output_file: Path,
    residue_output_file: Path,
    residue_bfactor_pdb: Path,
) -> None:
    wt_eval = read_eval_file(wt_eval_file)
    failed: list[str] = []
    mutation_impact_scores: dict[tuple[str, str, str], float] = {}
    max_scores_by_residue: dict[tuple[str, str], float] = {}

    for task in tasks:
        prefix = f"{task.chain}{task.resid}{task.to_resname}"
        mut_eval_file = find_or_collect_eval_output(output_dir, prefix)
        if not mut_eval_file.exists():
            failed.append(f"{prefix}: missing {mut_eval_file}")
            continue
        try:
            mut_eval = read_eval_file(mut_eval_file)
            signed_delta = vibrational_entropy_score(wt_eval, mut_eval, nmodes)
        except ValueError as exc:
            failed.append(f"{prefix}: {exc}")
            continue
        mutation_score = abs(signed_delta)
        mutation_impact_scores[(task.chain, task.resid, task.to_resname)] = mutation_score
        residue_key = (task.chain, task.resid)
        max_scores_by_residue[residue_key] = max(max_scores_by_residue.get(residue_key, 0.0), mutation_score)

    residue_raw_scores = [
        max_scores_by_residue.get((residue.chain, residue.resid), math.nan)
        for residue in residues
    ]
    valid_mask = residue_valid_mask(residues)
    smoothing_input = [0.0 if math.isnan(score) else score for score in residue_raw_scores]
    residue_smoothed_scores = gaussian_smooth_scores_self_include(
        residues,
        smoothing_input,
        valid_mask,
    )
    residue_score_map = {
        (residue.chain, residue.resid): score
        for residue, score in zip(residues, residue_smoothed_scores)
    }

    with mutation_output_file.open("w") as handle:
        handle.write("chain\tresid\twt_residue\tmut_residue\tmutation_impact_score(|ΔS|)\n")
        for task in tasks:
            mutation_score = mutation_impact_scores.get((task.chain, task.resid, task.to_resname))
            if mutation_score is None:
                continue
            handle.write(
                f"{task.chain}\t{task.resid}\t{task.from_resname}\t"
                f"{task.to_resname}\t{mutation_score:.8f}\n"
            )

    with residue_output_file.open("w") as handle:
        handle.write("chain\tresid\twt_residue\tresidue_impact_score(Ω)\n")
        for residue in residues:
            residue_score = residue_score_map.get((residue.chain, residue.resid), math.nan)
            residue_score_text = "NA" if math.isnan(residue_score) else f"{residue_score:.8f}"
            handle.write(
                f"{residue.chain}\t{residue.resid}\t{residue.resname}\t{residue_score_text}\n"
            )

    visualization_score_map = make_visualization_score_map(residue_score_map)
    write_residue_impact_bfactor_pdb(input_pdb, residue_bfactor_pdb, visualization_score_map)

    failed_file = mutation_output_file.with_suffix(".failed.txt")
    failed_file.write_text("\n".join(failed) + ("\n" if failed else ""))


def main() -> None:
    args = parse_args()
    pdb_id = args.input_pdb.stem
    output_dir = args.output_dir or Path(f"{pdb_id}_{args.chain}")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input PDB: {args.input_pdb}")
    print(f"Protein ID: {pdb_id}")
    print(f"Chain(s): {args.chain}")
    print(f"Output directory: {output_dir}")
    if args.entropy_mode is not None:
        print(f"Entropy modes: {args.entropy_mode}")
    else:
        print(f"Entropy mode percentage: {args.entropy_percent:g}%")
    print(f"Parallel jobs: {args.jobs}")
    print(f"ICed-ENM cores per NMA job: {args.nma_core}")
    print(f"Max mutation attempts: {args.max_attempts}")
    print(f"Base seed: {args.seed}")
    print(f"ICed-ENM executable: {args.iced_enm_bin}")
    if args.modeller_path is not None:
        print(f"Modeller path: {args.modeller_path}")
    else:
        print("Modeller path: active Python environment")
    print("Mutation generator: built-in Modeller routine")
    print(f"Amino acid scan set: {' '.join(AMINO_TYPES)}")

    wt_model = output_dir / f"{pdb_id}.pdb"
    shutil.copy2(args.input_pdb, wt_model)

    ref_ca = run_wt_iced_enm(args, output_dir)
    residues = read_ca_residues(ref_ca)
    nmode = determine_entropy_nmodes(args, len(residues))
    print(f"Normal modes for NMA: {nmode}")
    ref_ca = run_wt_iced_enm(args, output_dir, nmodes=nmode)
    cleanup_nma_outputs(output_dir, "wt")
    tasks = make_mutation_tasks(
        pdb_id=pdb_id,
        analysis_chains=args.chain,
        residues=residues,
        output_dir=output_dir,
        modeller_path=args.modeller_path.resolve() if args.modeller_path is not None else None,
        seed=args.seed,
        max_attempts=args.max_attempts,
    )
    print(f"Mutation residues: {len(residues)}")
    print(f"Mutation tasks: {len(tasks)}")
    log_output_dir = logs_dir(output_dir)
    generate_all_mutants(tasks, jobs=args.jobs, failed_list=log_output_dir / "failed_mutations.txt")
    run_all_mutant_nma(
        tasks,
        iced_enm_bin=args.iced_enm_bin.resolve(),
        nmodes=nmode,
        nma_core=args.nma_core,
        jobs=args.jobs,
        failed_list=log_output_dir / "failed_mutant_nma.txt",
    )
    write_mutation_impact_table(
        tasks,
        residues=residues,
        output_dir=output_dir,
        input_pdb=args.input_pdb,
        wt_eval_file=find_or_collect_eval_output(output_dir, "wt"),
        nmodes=nmode,
        mutation_output_file=output_dir / "mutation_impact_scores.tsv",
        residue_output_file=output_dir / "residue_impact_scores.tsv",
        residue_bfactor_pdb=output_dir / "residue_impact_bfactor.pdb",
    )


if __name__ == "__main__":
    main()
