import json

f = open("params.json")
params_dict = json.load(f)
f.close()

sensible_defaults = params_dict["sensible_defaults"]

# Input files
covariate_file = params_dict["covariate_file"]
genotype_file = params_dict["genotype_file"]
phenotype_file = params_dict["phenotype_file"]

population_name = params_dict["population_name"]
genotype_name = params_dict["genotype_name"]
covariate_name = params_dict["covariate_name"]
phenotype_name = params_dict["phenotype_name"]
model_type = params_dict["model_type"]

# Additional constants
chromosomes = params_dict["chromosomes"]
chunks = params_dict["chunks"]
window_size = params_dict["window_size"]
qsub_array = params_dict["qsub_array"]

def readable_bp_number(bp):
    if ((bp % 1e6) == 0):
        return f"{int(bp/1e6)}Mbase"
    elif ((bp % 1e3) == 0):
        return f"{int(bp/1e3)}kbase"
    else:
        return f"{bp}base"

model = f"{readable_bp_number(window_size)}_{model_type}"

if (sensible_defaults):
    out_covariate_dir = f"data/{population_name}/{genotype_name}_{phenotype_name}_{covariate_name}/covariates/"
    out_ped_dir = f"data/{population_name}/{genotype_name}/ped/"
    out_grm_dir = f"data/{population_name}/{genotype_name}/grm/"
    out_phenotype_dir = f"data/{population_name}/{phenotype_name}/"
    out_greml_intermediate_dir = f"data/{population_name}/{genotype_name}_{phenotype_name}_{covariate_name}/greml_intermediate/{model}/"
    out_hsq_dir = f"data/{population_name}/{genotype_name}_{phenotype_name}_{covariate_name}/hsq/{model}/"
    out_results_dir = f"data/{population_name}/{genotype_name}_{phenotype_name}_{covariate_name}/results/"
else: # Manually select output directories
    out_covariate_dir = "data/tri1_sQTL_covariates/"
    out_ped_dir = "data/tri1_sQTL_ped/"
    out_grm_dir = "data/tri1_sQTL_grm/"
    out_phenotype_dir = "data/tri1_sQTL_phenotype/"
    out_greml_intermediate_dir = "data/tri1_sQTL_greml_intermediate/"
    out_hsq_dir = "data/tri1_sQTL_hsq/"
    out_results_dir = "data/tri1_results/"

# Making output directories
import os

def makedirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

makedirs(out_covariate_dir)
makedirs(out_ped_dir)
makedirs(out_grm_dir)
makedirs(out_phenotype_dir)
makedirs(out_greml_intermediate_dir)
makedirs(out_hsq_dir)
makedirs(out_results_dir)

# Functions

# Rules
rule all:
    input:
        out_results_dir + model + ".csv"

# Processing genotype files
rule vcf_to_bed:
    input:
        genotype_file
    output:
        out_ped_dir + "genotype.bed",
        out_ped_dir + "genotype.bim",
        out_ped_dir + "genotype.fam",
        out_ped_dir + "genotype.log",
        out_ped_dir + "genotype.nosex"
    shell:
        "plink "
        "--vcf " + genotype_file + " "
        "--recode --const-fid --make-bed --mac 1 "
        "--out " + out_ped_dir + "genotype"

rule chr_grm:
    input:
        out_ped_dir + "genotype.bed",
        out_ped_dir + "genotype.bim",
        out_ped_dir + "genotype.fam",
        out_ped_dir + "genotype.log",
        out_ped_dir + "genotype.nosex"
    output:
        out_grm_dir + "genotype{chrom}.grm.N.bin",
        out_grm_dir + "genotype{chrom}.grm.bin",
        out_grm_dir + "genotype{chrom}.grm.id",
        out_grm_dir + "genotype{chrom}.log"
    shell:
        "gcta64 "
        "--bfile " + out_ped_dir + "genotype" + " "
        "--make-grm-bin --make-grm-alg 0 --chr {wildcards.chrom} "
        "--out " + out_grm_dir + "genotype{wildcards.chrom}"
 
rule whole_grm:
    input:
        expand(out_grm_dir + "genotype{chrom}.grm.N.bin", chrom=chromosomes),
        expand(out_grm_dir + "genotype{chrom}.grm.bin", chrom=chromosomes),
        expand(out_grm_dir + "genotype{chrom}.grm.id", chrom=chromosomes),
        expand(out_grm_dir + "genotype{chrom}.log", chrom=chromosomes)
    output:
        chr_list=out_grm_dir + "all_chrs_list.txt",
        grm_N=out_grm_dir + "complete.grm.N.bin",
        grm_bin=out_grm_dir + "complete.grm.bin",
        grm_id=out_grm_dir + "complete.grm.id",
        grm_log=out_grm_dir + "complete.log"
    shell:
        "echo > {output.chr_list} && " + \
        " ".join(["echo '" + out_grm_dir + "genotype" + str(chrom) + "' >> {output.chr_list} && " for chrom in chromosomes]) + \
        "gcta64 "
        "--make-grm-bin "
        "--mgrm {output.chr_list} "
        "--out " + out_grm_dir + "complete"

rule unrelated_grm:
    input:
        out_grm_dir + "complete.grm.N.bin",
        out_grm_dir + "complete.grm.bin",
        out_grm_dir + "complete.grm.id",
        out_grm_dir + "complete.log",
    output:
        out_grm_dir + "unrelated.grm.N.bin",
        out_grm_dir + "unrelated.grm.bin",
        out_grm_dir + "unrelated.grm.id",
        out_grm_dir + "unrelated.log"
    shell:
        "gcta64 "
        "--grm " + out_grm_dir + "complete "
        "--grm-cutoff 0.05 --make-grm-bin "
        "--out " + out_grm_dir + "unrelated"

rule excluded_samples_list:
    input:
        complete=out_grm_dir + "complete.grm.id",
        unrelated=out_grm_dir + "unrelated.grm.id",
    output:
        out_grm_dir + "excluded_samples.txt"
    shell:
        "diff "
        "--new-line-format=%L "
        "--unchanged-line-format='' "
        "{input.complete} "
        "{input.unrelated} | "
        r"sed 's/\t/ /g' > "
        "{output} || :"

# Processing phenotype files

rule bed_to_pheno:
    input:
        phenotype_file
    output:
        out_phenotype_dir + "phenotype.txt",
        out_phenotype_dir + "phenotype_info.txt"
    script:
        "scripts/bed_to_phenotype.jl"

checkpoint chunk_phenotype_info:
    input:
        out_phenotype_dir + "phenotype_info.txt"
    output:
        expand(out_phenotype_dir + "phenotype_info.txt_{chunk}", chunk = [str(chunk).rjust(7, '0') for chunk in [*range(0, chunks)]])
    shell:
        "split -n l/" + str(chunks) + " -a 7 -d {input} " + out_phenotype_dir + "phenotype_info.txt_"

# Splitting covariate files

rule split_covariate:
    input:
        covariate_file
    output:
        out_covariate_dir + "quant_cov.txt",
        out_covariate_dir + "qual_cov.txt"
    script:
        "scripts/fetal_scripts/process_covariates.jl"

# Checkpoint before running GREML

# Running GREML

if (model_type == "window"):
    greml_script = "scripts/run_window.jl"
    compile_script = "scripts/hsq_window.jl"
elif (model_type == "cis"):
    greml_script = "scripts/run_cis_greml_wrapper.jl"
    compile_script = "scripts/hsq_cis.jl"

if qsub_array:
    rule greml:
        input:
            genotypebed=out_ped_dir + "genotype.bed",
            genotypebim=out_ped_dir + "genotype.bim",
            genotypefam=out_ped_dir + "genotype.fam",
            genotypelog=out_ped_dir + "genotype.log",
            genotypenosex=out_ped_dir + "genotype.nosex",
            excludedsamples=out_grm_dir + "excluded_samples.txt",
            phenotype=out_phenotype_dir + "phenotype.txt",
            phenotypeinfo=expand(out_phenotype_dir + "phenotype_info.txt_{chunk}", chunk = [str(chunk).rjust(7, '0') for chunk in [*range(0, chunks)]]),
            quantcov=out_covariate_dir + "quant_cov.txt",
            qualcov=out_covariate_dir + "qual_cov.txt",
            grmN=expand(out_grm_dir + "genotype{chrom}.grm.N.bin", chrom = chromosomes),
            grmbin=expand(out_grm_dir + "genotype{chrom}.grm.bin", chrom = chromosomes),
            grmid=expand(out_grm_dir + "genotype{chrom}.grm.id", chrom = chromosomes),
            grmlog=expand(out_grm_dir + "genotype{chrom}.log", chrom = chromosomes)
        output:
            expand(out_greml_intermediate_dir + ".{chunk}_chunk.done", chunk = [str(chunk).rjust(7, '0') for chunk in [*range(0, chunks)]])
        params:
            genotypeprefix=out_ped_dir + "genotype",
            grmprefix=out_grm_dir + "genotype",
            outgremlintermediatedir=out_greml_intermediate_dir,
            outhsqdir=out_hsq_dir,
            windowsize = window_size
        run:
            with open(out_greml_intermediate_dir + "snakemake_info.json", "w") as outfile:
                json.dump({
                    "input": {
                        "excludedsamples": input.excludedsamples,
                        "phenotype": input.phenotype,
                        "phenotypeinfo": input.phenotypeinfo,
                        "quantcov": input.quantcov,
                        "qualcov": input.qualcov
                        },
                    "output": output,
                    "params": {
                        "genotypeprefix": params.genotypeprefix,
                        "grmprefix": params.grmprefix,
                        "outgremlintermediatedir": params.outgremlintermediatedir,
                        "outhsqdir": params.outhsqdir,
                        "windowsize": params.windowsize
                        }
                    },
                    outfile)
            shell("qsub "
            "-N qsub_greml "
            "-o logs/qsub_greml.log "
            "-e logs/qsub_greml.err "
            "-t 1:{chunks} "
            "scripts/run_cis_qsub.sh " + out_greml_intermediate_dir + "snakemake_info.json")
            shell("while [[ -n $(qstat | grep qsub_greml ) ]]; do "
                "sleep 300; "
                "done")
            shell("rm " + out_greml_intermediate_dir + "snakemake_info.json")
else:
    rule greml:
        input:
            genotypebed=out_ped_dir + "genotype.bed",
            genotypebim=out_ped_dir + "genotype.bim",
            genotypefam=out_ped_dir + "genotype.fam",
            genotypelog=out_ped_dir + "genotype.log",
            genotypenosex=out_ped_dir + "genotype.nosex",
            excludedsamples=out_grm_dir + "excluded_samples.txt",
            phenotype=out_phenotype_dir + "phenotype.txt",
            phenotypeinfo=out_phenotype_dir + "phenotype_info.txt_{chunk}",
            quantcov=out_covariate_dir + "quant_cov.txt",
            qualcov=out_covariate_dir + "qual_cov.txt",
            grmN=expand(out_grm_dir + "genotype{chrom}.grm.N.bin", chrom = chromosomes),
            grmbin=expand(out_grm_dir + "genotype{chrom}.grm.bin", chrom = chromosomes),
            grmid=expand(out_grm_dir + "genotype{chrom}.grm.id", chrom = chromosomes),
            grmlog=expand(out_grm_dir + "genotype{chrom}.log", chrom = chromosomes)
        output:
            out_greml_intermediate_dir + ".{chunk}_chunk.done"
        params:
            genotypeprefix=out_ped_dir + "genotype",
            grmprefix=out_grm_dir + "genotype",
            outgremlintermediatedir=out_greml_intermediate_dir,
            outhsqdir=out_hsq_dir,
            windowsize = window_size
        script:
            greml_script

rule compile_greml:
    input:
        expand(out_greml_intermediate_dir + ".{chunk}_chunk.done", chunk = [str(chunk).rjust(7, '0') for chunk in [*range(0, chunks)]])
    output:
        out_results_dir + model + ".csv"
    params:
        outhsqdir=out_hsq_dir
    script:
        compile_script

