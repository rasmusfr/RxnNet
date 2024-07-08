<h1 align="center">
  <picture>
    <img alt="Logo" src="resources/logo.png"
height="70">
  </picture>
</h1>

<h4 align="center">

[![Requires Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)

</h4>

RxnNet is a reaction network based python package for property prediction in solids.

## Installation
Install from source:

```sh
git clone https://github.com/rasmusfr/rxnnet.git
cd rxnnet
pip install .
```

## Usage

### Reaction network generation
Using the default database, RxnNet can generate a reaction network (i.e. a series of balanced reactions). 
At the moment, reactions with up to 4 compounds are supported.

```python
from rxnnet.generate_reactions import GenerateRN

target_id = 'mp-23193'  # target compound materials project id (mp-23193 is KCl)

reaction_network_gen = GenerateRN(target_id)
rn_data = reaction_network_gen.balanced_reactions(save=True)

print(rn_data[:50].to_markdown())
```
### Predictions from the reaction network
From the generated reaction network we can obtain a prediction for the formation enthalpy.

```python
import pandas as pd

from rxnnet.evaluate_reactions import EvaluateRN

target_id = 'mp-23193'

preferred_method = 'e_r2scan'  # preferred method; r^2SCAN electronic energy
fallback_method = 'e_gga_gga_u'  # fallback method; GGA/GGA+U electronic energy
reference_method = 'hf_ref'  # reference method; NBS enthalpy of formation
mode = 'ssw+cf'  # calculation mode; structural similarity weighting + chemistry filter

reaction_network_eval = EvaluateRN(target_id=target_id, preferred_method=preferred_method, reference_method=reference_method,
                                   rn=rf'user_reactions/{target_id}.pkl.gz', mode=mode, fallback_method=fallback_method)

results = reaction_network_eval.rn_evaluate()
print(pd.DataFrame(results)[:50].to_markdown())
```

## Example notebooks
| Notebook                               | Description                                        |
|----------------------------------------|----------------------------------------------------|
| [Intro](examples/intro-mp-23193.ipynb) | Introduction to reaction generation and prediction |
## Datasets
Datasets and results related to RxnNet can be found [here](https://doi.org/10.11583/DTU.25897420).
## License
This repository is licensed under the [MIT license](LICENSE)
## Citing RxnNet
```bib
This repository was created by Rasmus Fromsejer (Technical University of Denmark) to supplement the research paper "Accurate Formation Enthalpies of Solids Using Reaction Networks" in npj computational materials by Rasmus Fromsejer, Bjørn Maribo-Mogensen, Georgios Kontogeorgis and Xiaodong Liang (accepted in principle).
```
## Acknowledgement
The author wishes to thank the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program (Grant Agreement no. 832460), ERC Advanced Grant project “New Paradigm in Electrolyte Thermodynamics” and the Department of Chemical Engineering at the Technical University of Denmark for funding this research.
