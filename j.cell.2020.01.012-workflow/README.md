# Workflow to replicate [Cell Data](https://doi.org/10.1016/j.cell.2020.01.012)

To execute, run `snakemake --profile local_run`.  This will perform the
analysis pipeline with 2 wget isntances.  Final calls will be produced in the
`ibd_final` directory under the specified `working_directory`.
