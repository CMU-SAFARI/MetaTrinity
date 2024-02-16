# MetaFast: Enabling Fast Metagenomic Classification via Seed Counting and Edit Distance Approximation

We introduce MetaFast, a tool based on heuristics, to achieve a fundamental improvement in accuracy-runtime tradeoff over existing methods. MetaFast delivers accuracy comparable to the alignment-based and highly accurate tool Metalign but with significantly enhanced efficiency. In MetaFast, we accelerate memory-frugal reference database indexing and filtering. We further employ heuristics to accelerate read mapping. Our evaluation demonstrates that MetaFast achieves a 3.9x speedup over Metalign without compromising accuracy.  
MetaFast has been rigorously debugged and validated using GCC 12.1.0-16.

Described by Arvid E. Gollwitzer et al. (current version at https://doi.org/10.48550/arXiv.2311.02029).

## Table of Contents
- [Setup the Database and Index Structure ](#install)
- [Getting Help](#contact)
- [Citing MetaFast](#cite)


## <a name="install"></a>Installation
```sh
# TBD: Instructions to get fna files and make MMIs!
```


##  <a name="contact"></a>Getting Help
If you have suggestions for improvement, new applications, or wish to collaborate, please contact arvid.gollwitzer@safari.ethz.ch.  
If you encounter bugs or have further questions or requests, you may raise an issue on the issue page.


## <a name="cite"></a>Citing MetaFast
If you use MetaFast in your work, please cite:

> Gollwitzer, Arvid E., et al. "MetaFast: Enabling Fast Metagenomic Classification via Seed Counting and Edit Distance Approximation." 
> arXiv:2311.02029 (2023).

You may use the following BibTeX:

```bibtex
@misc{gollwitzer2023metafast,
      title={MetaFast: Enabling Fast Metagenomic Classification via Seed Counting and Edit Distance Approximation}, 
      author={Arvid E. Gollwitzer and Mohammed Alser and Joel Bergtholdt and Joel Lindegger and Maximilian-David Rumpf and Can Firtina and Serghei Mangul and Onur Mutlu},
      year={2023},
      eprint={2311.02029},
      archivePrefix={arXiv},
      primaryClass={q-bio.GN}
}

