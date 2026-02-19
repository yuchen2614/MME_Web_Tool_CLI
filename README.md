# Perfect_match & Blastp Analysis CLI

Command-line pipeline for comparing input query proteins and the human reference proteome

## Installation

```bash
git clone https://github.com/yuchen2614/MME_Web_Tool_CLI.git
```

## Environment setting

```bash
cd MME_Web_Tool_CLI
pip install -r requirements.txt
```

## Usage

1. Edit file: MME_Web_Tool_CLI/analysis_pipeline_code/CLI_config.yaml
(choose a "Matching_Logic" then change corresponding query_fasta_path and its parameter)
2. Make terminal location to see "MME_Web_Tool_CLI" folder.
3. Then use terminal command "python -m MME_Web_Tool_CLI.analysis_pipeline_code.CLI_run" to run the analysis.

```bash
python -m MME_Web_Tool_CLI.analysis_pipeline_code.CLI_run
```
