[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6308926.svg)](https://doi.org/10.5281/zenodo.6308926)

# julia-recon-advances-in-spiral-fmri
Example pipeline for undersampled (R=4) spiral image reconstruction with static B0 correction 
- using [MRIReco.jl](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/), an MRI reconstruction framework written in Julia ([Knopp and Grosser, MRM, 2021](https://doi.org/10.1002/mrm.28792))
- and raw data (ISMRMRD, NIfTI) publicly available for download within the [ETH Research Collection](https://doi.org/10.3929/ethz-b-000487412) and described as a [Data in Brief](https://doi.org/10.1016/j.dib.2022.108050) article.

This repository showcases an iterative image reconstruction (conjugate gradient SENSE) of undersampled high-resolution single-shot spiral fMRI data, including accelerated static B0 correction using pre-acquired and -processed B0 and sensitivity maps. The dataset is described in detail in a [Data in Brief](https://doi.org/10.1016/j.dib.2022.108050) article. The spiral imaging results from that data, reconstructed with a similar algorithm (in Matlab), were published in our paper [Advances in spiral fMRI: A high-resolution study with single-shot acquisition](https://doi.org/10.1016/j.neuroimage.2021.118738) (Kasper et al., NeuroImage, 2022).


# Related Work

## Publications

- Data: [Advances in spiral fMRI: A high-resolution dataset](https://doi.org/10.1016/j.dib.2022.108050), Data In Brief, 2022
- Results: [Advances in spiral fMRI: A high-resolution study with single-shot acquisition](https://doi.org/10.1016/j.neuroimage.2021.118738), NeuroImage, 2022

## Code

- [Analysis and Representation Code (Matlab)](https://github.com/mrikasper/paper-advances-in-spiral-fmri)

## Data

- [Data Collection for the Article "Advances in Spiral fMRI - A High-resolution Study with Single-shot Acquisition"](https://doi.org/10.3929/ethz-b-000487412), ETH Research collection, 2021
    - Additional Metadata and data description: Article [Advances in spiral fMRI: A high-resolution dataset](https://doi.org/10.1016/j.dib.2022.108050), Data in Brief, 2022
