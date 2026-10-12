#!/bin/bash

# Submit reproducible DREAMSim optical-photon jobs on TTU HPCC.
# Adapted from Samuel McKinley's dSiPM_SPAD_NN/batch_Sims.sh.

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
campaign_file=${1:-${repo_dir}/production/campaigns/cpad2026_pilot_100gev.tsv}

required_dreamsim_commit=54b2447c323dd7eb5093fda5baf9ee433bad1bf6
dreamsim_dir=${DREAMSIM_DIR:-${repo_dir}/../DREAMSim}
sim_build_dir=${SIM_BUILD_DIR:-${dreamsim_dir}/sim/build}
sim_executable=${SIM_EXECUTABLE:-${sim_build_dir}/exampleB4b}
sim_macro=${SIM_MACRO:-${dreamsim_dir}/sim/paramBatch03_single.mac}
singularity_image=${SINGULARITY_IMAGE:-/lustre/research/hep/yofeng/SimulationEnv/alma9forgeant4_sbox}
output_base=${OUTPUT_BASE:-/lustre/work/${USER}/flash_dreamsim/cpad2026_pilot_100gev}
partition=${PARTITION:-nocona}
memory=${MEMORY:-16G}
time_limit=${TIME_LIMIT:-12:00:00}
seed_base=${SEED_BASE:-202610110}
submit=${SUBMIT:-0}
validator_python=${VALIDATOR_PYTHON:-python3}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

if [ ! -f "$campaign_file" ]; then
    echo "ERROR: Missing campaign file: $campaign_file" >&2
    exit 2
fi
if [ ! -d "$dreamsim_dir/.git" ]; then
    echo "ERROR: DREAMSim is not a git checkout: $dreamsim_dir" >&2
    exit 2
fi

if [ ! -d "$repo_dir/.git" ]; then
    echo "ERROR: dSiPM is not a git checkout: $repo_dir" >&2
    exit 2
fi

dsipm_commit=$(git -C "$repo_dir" rev-parse HEAD)
if ! git -C "$repo_dir" diff --quiet || ! git -C "$repo_dir" diff --cached --quiet; then
    echo "ERROR: dSiPM has tracked changes. Submit from a committed checkout." >&2
    exit 2
fi

actual_commit=$(git -C "$dreamsim_dir" rev-parse HEAD)
if [ "$actual_commit" != "$required_dreamsim_commit" ]; then
    echo "ERROR: DREAMSim must be at $required_dreamsim_commit" >&2
    echo "       Found $actual_commit in $dreamsim_dir" >&2
    exit 2
fi
if ! git -C "$dreamsim_dir" diff --quiet || ! git -C "$dreamsim_dir" diff --cached --quiet; then
    echo "ERROR: DREAMSim has tracked source changes. Use a clean pinned checkout." >&2
    exit 2
fi
if [ ! -s "$sim_macro" ]; then
    echo "ERROR: Missing simulation macro: $sim_macro" >&2
    exit 2
fi

if [ "$submit" = "1" ]; then
    command -v sbatch >/dev/null 2>&1 || {
        echo "ERROR: sbatch is unavailable; run this on TTU HPCC." >&2
        exit 2
    }
    if [ ! -x "$sim_executable" ]; then
        echo "ERROR: Missing DREAMSim executable: $sim_executable" >&2
        exit 2
    fi
    if [ ! -e "$singularity_image" ]; then
        echo "ERROR: Missing Singularity image: $singularity_image" >&2
        exit 2
    fi
    if ! "$validator_python" -c 'import ROOT' >/dev/null 2>&1 \
        && ! "$validator_python" -c 'import uproot, awkward' >/dev/null 2>&1; then
        echo "ERROR: $validator_python has neither PyROOT nor uproot/awkward." >&2
        exit 2
    fi
fi

launcher_sha256=$(sha256_file "$repo_dir/production/submit_dreamsim_campaign.sh")
validator_sha256=$(sha256_file "$repo_dir/production/validate_dreamsim_output.py")
campaign_sha256=$(sha256_file "$campaign_file")
sim_executable_sha256="not-checked-dry-run"
if [ -x "$sim_executable" ]; then
    sim_executable_sha256=$(sha256_file "$sim_executable")
fi

mkdir -p "$output_base/jobs" "$output_base/logs" "$output_base/provenance"
for source_file in \
    "$repo_dir/production/submit_dreamsim_campaign.sh" \
    "$repo_dir/production/validate_dreamsim_output.py" \
    "$campaign_file"; do
    snapshot="$output_base/provenance/$(basename "$source_file")"
    if [ -e "$snapshot" ] && ! cmp -s "$source_file" "$snapshot"; then
        echo "ERROR: Provenance snapshot differs from current source: $snapshot" >&2
        echo "       Use a new OUTPUT_BASE for a changed campaign." >&2
        exit 2
    fi
    cp "$source_file" "$snapshot"
done

manifest="$output_base/submission_manifest.tsv"
if [ ! -e "$manifest" ]; then
    printf 'sample\tparticle\tenergy_gev\tjob_index\tevents\tseed_1\tseed_2\tslurm_job_id\tdreamsim_commit\tdsipm_commit\tmacro_sha256\tcampaign_sha256\tlauncher_sha256\tvalidator_sha256\tsim_executable_sha256\n' > "$manifest"
fi

row_index=0
while IFS=$'\t' read -r sample particle energy_gev jobs events_per_job extra; do
    if [ -z "${sample:-}" ] || [[ "$sample" == \#* ]]; then
        continue
    fi
    if [ -n "${extra:-}" ]; then
        echo "ERROR: Too many columns in campaign row for $sample" >&2
        exit 2
    fi
    if ! [[ "$jobs" =~ ^[1-9][0-9]*$ && "$events_per_job" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: jobs and events_per_job must be positive integers for $sample" >&2
        exit 2
    fi
    row_index=$((row_index + 1))
    sample_dir="$output_base/$sample"
    mkdir -p "$sample_dir"

    for ((job_index=0; job_index<jobs; job_index++)); do
        job_tag=$(printf '%03d' "$job_index")
        seed_1=$((seed_base + row_index * 100000 + job_index * 2))
        seed_2=$((seed_1 + 1))
        job_name="FLASH_${sample}_${job_tag}"
        run_number=$row_index
        run_seq=$job_index
        expected_root="${sample_dir}/mc_${job_name}_run${run_number}_${run_seq}_FullDetector_${events_per_job}evt_${particle}_${energy_gev}_${energy_gev}.root"
        job_script="$output_base/jobs/${job_name}.sh"
        job_macro="$output_base/jobs/${job_name}.mac"
        job_validator="$output_base/jobs/${job_name}.validator.py"

        if [ -s "$expected_root.validation.json" ]; then
            echo "SKIP: validated output already exists: $expected_root"
            continue
        fi

        {
            echo "/random/setSeeds ${seed_1} ${seed_2}"
            awk -v root_pre="${sample_dir}/mc" '
                /^[[:space:]]*#\$\$\$[[:space:]]+rootPre[[:space:]]/ {
                    print "#$$$ rootPre       " root_pre
                    next
                }
                /^[[:space:]]*\/physics_list\/list[[:space:]]*$/ { next }
                { print }
            ' "$sim_macro"
        } > "$job_macro"
        job_macro_sha256=$(sha256_file "$job_macro")
        cp "$repo_dir/production/validate_dreamsim_output.py" "$job_validator"

        cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH -J ${job_name}
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p ${partition}
#SBATCH --mem=${memory}
#SBATCH --time=${time_limit}
#SBATCH -o ${output_base}/logs/%x.%j.out
#SBATCH -e ${output_base}/logs/%x.%j.err

set -euo pipefail

output_dir=${sample_dir}
tmp_dir=/lustre/scratch/\${USER}/flash_\${SLURM_JOB_ID}
mkdir -p "\$output_dir" "\$tmp_dir/home" "\$tmp_dir/xdg-cache" "\$tmp_dir/singularity-cache" "\$tmp_dir/singularity-tmp"
export TMPDIR="\$tmp_dir"
export SINGULARITY_CACHEDIR="\$tmp_dir/singularity-cache"
export SINGULARITY_TMPDIR="\$tmp_dir/singularity-tmp"
trap 'rm -rf "\$tmp_dir"' EXIT
ulimit -c 0 || true

macro="\$tmp_dir/${job_name}.mac"
cp ${job_macro} "\$macro"

cd "\$output_dir"
singularity exec --cleanenv --bind /lustre:/lustre --bind "\$tmp_dir":/tmp ${singularity_image} bash --noprofile --norc -c 'export HOME=/tmp/home XDG_CACHE_HOME=/tmp/xdg-cache ROOT_HIST=0; ${sim_executable} -b '"\$macro"' -numberOfEvents ${events_per_job} -eventsInNtupe ${events_per_job} -jobName ${job_name} -runNumber ${run_number} -runSeq ${run_seq} -runConfig FullDetector -gun_particle ${particle} -gun_energy_min ${energy_gev} -gun_energy_max ${energy_gev} -createNtuple true -miniNtuple false -saveTruthHits true -sipmType 1'

root_file=${expected_root}
if [ ! -s "\$root_file" ]; then
    echo "ERROR: DREAMSim did not produce the expected ROOT file: \$root_file" >&2
    exit 2
fi

${validator_python} ${job_validator} "\$root_file" --expected-events ${events_per_job} --expected-particle '${particle}' --expected-energy-gev ${energy_gev} --dreamsim-commit ${required_dreamsim_commit} --dsipm-commit ${dsipm_commit} --macro-sha256 ${job_macro_sha256} --campaign-sha256 ${campaign_sha256} --launcher-sha256 ${launcher_sha256} --validator-sha256 ${validator_sha256} --sim-executable-sha256 ${sim_executable_sha256} --seed-1 ${seed_1} --seed-2 ${seed_2}
EOF
        chmod +x "$job_script"

        slurm_job_id="DRY_RUN"
        if [ "$submit" = "1" ]; then
            slurm_job_id=$(sbatch --parsable "$job_script")
            echo "Submitted $job_name as SLURM job $slurm_job_id"
        else
            echo "Prepared $job_script"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$particle" "$energy_gev" "$job_index" "$events_per_job" "$seed_1" "$seed_2" "$slurm_job_id" "$required_dreamsim_commit" "$dsipm_commit" "$job_macro_sha256" "$campaign_sha256" "$launcher_sha256" "$validator_sha256" "$sim_executable_sha256" >> "$manifest"
    done
done < "$campaign_file"

echo "Campaign preparation complete. Manifest: $manifest"
if [ "$submit" != "1" ]; then
    echo "Dry run only. Set SUBMIT=1 on TTU HPCC to submit the prepared campaign."
fi
