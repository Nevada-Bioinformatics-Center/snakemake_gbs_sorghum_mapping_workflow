with open(snakemake.output[0], 'w') as writer:
    for sample in sorted(snakemake.params.samples):
        writer.write(sample + "\n")
