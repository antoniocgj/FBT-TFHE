# Revisiting the functional bootstrap in TFHE

This repository holds the source code for reproducing the results of [Revisiting the functional bootstrap in TFHE](https://doi.org/10.46586/tches.v2021.i2.229-253). See [instructions.md](https://github.com/antoniocgj/FBT-TFHE/blob/main/instructions.md) for building instructions.

## Paper Abstract

The FHEW cryptosystem introduced the idea that an arbitrary function can be evaluated within the bootstrap procedure as a table lookup. The faster bootstraps of TFHE strengthened this approach, which was later named Functional Bootstrap (Boura et al., CSCML’19). From then on, little effort has been made towards defining efficient ways of using it to implement functions with high precision. In this paper, we introduce two methods to combine multiple functional bootstraps to accelerate the evaluation of reasonably large look-up tables and highly precise functions. We thoroughly analyze and experimentally validate the error propagation in both methods, as well as in the functional bootstrap itself. We leverage the multi-value bootstrap of Carpov et al. (CT-RSA’19) to accelerate (single) lookup table evaluation, and we improve it by lowering the complexity of its error variance growth from quadratic to linear in the value of the output base. Compared to previous literature using TFHE’s functional bootstrap, our methods are up to 2.49 times faster than the lookup table evaluation of Carpov et al. (CT-RSA’19) and up to 3.19 times faster than the 32-bit integer comparison of Bourse et al. (CT-RSA’20). Compared to works using logic gates, we achieved speedups of up to 6.98, 8.74, and 3.55 times over 8-bit implementations of the functions ReLU, Addition, and Maximum, respectively.

## Citation

```json
@article{Guimarães_Borin_Aranha_2021, 
  title={Revisiting the functional bootstrap in TFHE},
  volume={2021}, 
  url={https://tches.iacr.org/index.php/TCHES/article/view/8793}, 
  DOI={10.46586/tches.v2021.i2.229-253}, 
  number={2}, 
  journal={IACR Transactions on Cryptographic Hardware and Embedded Systems}, 
  author={Guimarães, Antonio and Borin, Edson and Aranha, Diego F.}, 
  year={2021}, 
  month={Feb.}, 
  pages={229-253} 
}
```