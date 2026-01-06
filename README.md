### Conda enviroment:
This pipeline is compatible with the conda enviroment used in our [in-house clustering pipeline]([url](https://github.com/vijaybioinfo/clustering)), if its already available, please use it, otherwise follow these instructions:

```
conda create --name clustering python=3.8
conda activate clustering
pip install -r requirements.txt
conda deactivate
```

### Run pipeline:
Add the necessary information to the YAML file. Take the example_config.yaml file in the example folder as a reference.

After you've added the necessary information to the YAML file you can call the pipeline as:
sh Xenium-clustering/run.sh -y /path/to/your/yaml/example_config.yaml



