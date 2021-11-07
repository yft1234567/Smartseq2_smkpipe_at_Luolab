# import glob
# import os


# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")

# samples = (
#     pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
#     .set_index("sample_name", drop=False)
#     .sort_index()
# )

# validate(samples, schema="../schemas/samples.schema.yaml")

# def get_final_outputs():
#     # Counting from single batch, omitting down-sampling and saturation test
#     final_outputs = expand(
#         "workflow/data/{sample.user}/{sample.library}/outs/{sample.sample_name}_final_count",
#         sample = samples.itertuples()
#     )
#     return final_outputs

# def get_final_outputs():
#     if config["saturation_test"]["activate"]:
#         # Perform down-sampling and count on various levels of sequencing depths (units)
#         final_outputs = expand(
#             "workflow/data/{sample.user}/{sample.directory}/outs/{sample.sample_name}_{units}_final_counts",
#             sample = samples.itertuples(),
#             units = get_units()
#         )
#     else:
#         # Counting from single batch, omitting down-sampling and saturation test
#         final_outputs = expand(
#             "workflow/data/{sample.user}/{sample.directory}/outs/{sample.sample_name}_final_counts",
#             sample = sample.itertuples()
#         )
#     return final_outputs.to_list()


# def get_units():
#     max_unit = config['saturation_test']['max_unit']
#     step = config['saturation_test']['step']
#     units = list(range(step, max_unit, step))
#     units.append(max_unit)
#     return units

