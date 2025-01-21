# ISI_Wiretap_Channels.Code_Design (IN PREPARATION)

#### Design and Analysis of a Concatenated Code for Intersymbol Interference Wiretap Channels
> + Preprint is available at [arXiv:*2501.07561 [cs.IT]*](https://arxiv.org/abs/2501.07561).
> + See our [IEEE ISIT 2022](https://ieeexplore.ieee.org/abstract/document/9834578#citations) paper: Matched Information Rate Codes for Binary-Input Intersymbol Interference Wiretap Channels


#### Simulation Files:
> + Monte-Carlo (*Encoding/Decoding*) [Code Simulation](https://github.com/arianouri/ISI_Wiretap_Channels.Code_Design/tree/main/%5BSIMULATION_FILES%5D%20Code%20Design/MC_Encoding_Decoding)
> 1. `main_step1_random_PCM.m` Given degree distributions, this file generates parity-check matrices (PCMs) of the outer LDPC code stage without girth 4. (You need to have CVX installed + licensed MOSEK Optimization Tool version 10.)
> 2. `main_step2_uppertri_PCM.m` Given a PCM, this file produces the corresponding upper triangular format proper for reduced-complexity encoding of LDPC codes (according to Appendix A of *``Modern Coding Theory,''*  Book by Richardson and Urbanke).
> 3. `main_step3_MC_enc_dec` Monte-Carlo simulation of the encoder, the ISI wiretap channel, and the legitimate receiver's message passing decoder.
> + Inner Trellis Code Stage Design: *Files are in prep.*
> + Outer LDPC Code Stage Design: *Files are in prep.*
