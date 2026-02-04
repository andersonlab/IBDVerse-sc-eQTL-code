

params {
    RUN = "$RUN_ID"
    PATH2 = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/CELLRANGER_V7/rectum"
    // PATH2 = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/CELLRANGER_V7/TI"
    // PATH2 = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/CELLRANGER_V7/blood"
    input = 'cellbender' 
    lisi{
        run_process=false
    }
    replace_genotype_ids=false
    webtransfer = true
    write_h5=false
    encrypt = false
    rsync_to_web_file='/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/SETUP_fech_input_prep/scripts/rsync_to_web.sh'
    cohorts_to_drop_from_GT_Relatednes_check='GT_UKBB'
    // skip_qc=true
    tmpdir = "${launchDir}/tmp"
    project_name = 'Anderson_Lab'
    cellbender_location='/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/CELLRANGER_V7/rectum/cellbender/' //!!!!! if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'
    cellbender_resolution_to_use='0pt1'
    extra_metadata = "$PATH2/Extra_Metadata.tsv"
    extra_sample_metadata =""
    genotype_phenotype_mapping_file = ''
        // output_dir = outdir= "${launchDir}/results_rsync/results"
    output_dir = outdir= "${launchDir}/results"
    run_celltype_assignment=true
    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full tranche
    input_data_table = "$PATH2/input.tsv" //this has to be a full path
    // input_data_table = "$outdir/handover/Summary_plots/$RUN_ID/Fetch Pipeline/Input/input_table.tsv"
    // existing_cellsnp=''
    
    celltype_assignment{
        run_celltype_assignment=true
        run_azimuth=true
        run_keras=true
        run_celltypist=true
        run_scpred=false
    }

    celltype_prediction {
        keras {
            //# https://yascp.cog.sanger.ac.uk/public/singularity_images/eqtl_26_10_2022.img
            keras_model = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/data/keras_ti_model.h5'
            keras_weights_df = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/data/keras_ti_weights.tsv.gz'
            h5_layer = 'log1p_cp10k'
            keras_model_cluster_labels = "${projectDir}/assets/keras_cluster/data-cluster_labels.csv"
            filter_top_cell_probabilities = '0.5,0.75'
            save_all_probabilities = '--save_all_probabilities'
	    }
    }

    celltypist {
        run = true
        description = """https://github.com/Teichlab/celltypist"""
        remove_workdir = false
        sample_plot_probs = false
        copy_mode = "rellink"
        models = ['Immune_All_High.pkl','Immune_All_Low.pkl','Cells_Intestinal_Tract.pkl']
        //# OPTIONS:
        //# Adult_Mouse_Gut.pkl         COVID19_Immune_Landscape.pkl  Human_IPF_Lung.pkl    Lethal_COVID19_Lung.pkl
        //# Autopsy_COVID19_Lung.pkl    Developing_Human_Brain.pkl    Human_Lung_Atlas.pkl  models.json
        //# Cells_Fetal_Lung.pkl        Developing_Human_Thymus.pkl   Human_PF_Lung.pkl     Nuclei_Lung_Airway.pkl
        //# Cells_Intestinal_Tract.pkl  Developing_Mouse_Brain.pkl    Immune_All_High.pkl   Pan_Fetal_Human.pkl
        //# Cells_Lung_Airway.pkl       Healthy_COVID19_PBMC.pkl      Immune_All_Low.pkl
    }

	genotype_input {
        run_with_genotype_input=false
        vireo_with_gt=false
        posterior_assignment = false //if this is set to true, we will perform the genotype donor matching after the deconvolution is performed.
        subset_genotypes = false //if activated this will use th IDs provided to estimate the donors, otherwise it will match against full cohort
        full_vcf_file = '' //this could be a list of vcfs, in which case have to merge them 
        tsv_donor_panel_vcfs = ""
    }
    souporcell {
        run = false
    }

    harmony {
        run = false
    }
    filter_multiplets{
        run_process = true
        doubletDetection{
             run_process = false           
        }
        doubletDecon{
            run_process = false
        }
        scDblFinder{
            run_process = false
        }
        scds{
            run_process = false
        }
        doubletFinder{
            run_process = false
        }
        scrublet{
            run_process= true
        }
    }

    skip_preprocessing = false
    do_deconvolution = false

    harmony{
        run_process= false
        description = 'Parameters for harmony'
        variables_and_thetas{
            description = 'Tuples of metadata columns and corresponding thetas'
            value = [
                [ variable: 'experiment_id', theta: 1.0 ],
            ]
        }
    }

    bbknn{
        run_process = false
        description = 'Parameters for BBKNN'
        batch_variable{
            description = 'Variable to use for batch correction'
            value = 'experiment_id'
        }
    }

    lisi{
        run_process = false
        description = 'Parameters for Local Inverse Simpsons Index (LISI).'
        variables{
            description = 'Metadata variables to compute LISI over.'
            value = 'experiment_id'
        }
    }
}


process {

    withName: plot_distributions{
        containerOptions = "--containall --cleanenv --workdir /tmp -B /tmp"
    }

    withName: cellex_cluster_markers{
        maxForks=7
        memory = 300.GB
    }
    
    withName: GATHER_DATA{
        maxForks=7
        memory = 250.GB
    }
    withName: LISI{
        maxForks=7
        memory = 500.GB
    }
    withName: cluster_validate_resolution_keras{
        memory = 300.GB
    }

    withName: umap_calculate_and_plot{
        memory = 300.GB
    }

    withName: SPLIT_BATCH_H5AD{
        memory = { 100.GB * task.attempt }
    }

    withName: merge_samples_from_h5ad{
        memory = { 100.GB * task.attempt }
    }


    withName: sccaf_assess_clustering{
        memory = 300.GB
    }
    
}