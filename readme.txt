This program preprocess folders with genomes and run the EloE program. For each sample the EloE generates files where the proteins are ranked according to the level of translation elongation efficiency. Then this script devides those lists into a given number of quantiles, outputing the potentially highly expressed genes.

Before starting the program, you need to enter the path to EloE in the path.ini file according to the sample.
You can also add: input_path - a directory with folders for each sample containing .gbk files (one file per sample).
  output_path - empty directory (cleared before starting the program), where folders are added for each sample, containing the splitting of genes by quantiles.

You can set the input_path with the -i(--indir) command line option and the output_path with -o (--outdir) (mandatory if they weren't added to the path.ini file).
If any path is specified in both the command line and path.ini, only the path from the command line is used.
You can also set the number of quantiles with the -q(--quantiles) command, the default value is 4.
