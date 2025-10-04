# ISI_Wiretap_Channels.Code_Design

### Design and Analysis of a Concatenated Code for Intersymbol Interference Wiretap Channels

- Preprint available at [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561).  
- Presented at [IEEE ISIT 2022](https://ieeexplore.ieee.org/abstract/document/9834578).  

---
## Simulation Files:

### (I) Monte-Carlo: *Encoding/Decoding* [\[MC_Encoding_Decoding\]](https://github.com/arianouri/ISI_Wiretap_Channels.Code_Design/tree/main/%5BSIMULATION_FILES%5D%20Code%20Design/MC_Encoding_Decoding)

1. `\MC_Encoding_Decoding\main_step1_random_PCM.m`  
   * Given degree distributions, this file generates parity-check matrices (PCMs) for the outer LDPC code stage without girth-4.  
   * **Note:** Requires CVX and a licensed version of the MOSEK Optimization Tool (version 10).  

2. `\MC_Encoding_Decoding\main_step2_uppertri_PCM.m`  
   * Given a PCM, this file produces the corresponding upper-triangular format suitable for reduced-complexity encoding of LDPC codes.  
   * **Reference:** Appendix A of *Modern Coding Theory* by Richardson and Urbanke.  

3. `\MC_Encoding_Decoding\main_step3_MC_enc_dec.m`  
   * Performs a Monte-Carlo simulation of the encoder, the ISI wiretap channel, and the legitimate receiver’s message-passing decoder.  

---

### (II) Design: *Inner-Stage Trellis Code*  [\[S0--DESIGN-Inner Trellis Code Stage\]](https://github.com/arianouri/ISI_Wiretap_Channels.Code_Design/tree/b20cbda79aa2a426671f5203387c611d9c9f0814/%5BSIMULATION_FILES%5D%20Code%20Design/S0--DESIGN-Inner%20Trellis%20Code%20Stage)

> **Note:** You must have CVX installed along with a licensed version of the MOSEK Optimization Tool (version 10).

1. `\S0--DESIGN-Inner Trellis Code Stage\RULE_1_Expectation-Maximization\Main_Equi_no_MC_nth_45_30.m`  
   * **Rule 1:** This file optimizes a finite-order Markov source at the input of the ISI-WTC, modeling the received signals from a phased array.  
   * **Check** the following in [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561):  
     * Appendix A  
     * Examples 1 and 2  

2. `\S0--DESIGN-Inner Trellis Code Stage\RULE_2\Rule_2.m`  
   * **Rule 2:** This file solves the optimization problem stated in Eq. (12) of [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561).  

3. `\S0--DESIGN-Inner Trellis Code Stage\RULE_3\Rule_3.m`  
   * **Rule 3:** This file solves the optimization problem stated in Eq. (13) of [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561).  
   * **Check** the following in [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561):  
     * Criterion 1  
     * Example 3  

4. `\S0--DESIGN-Inner Trellis Code Stage\RULE_4\Rule_4.m`  
   * **Rule 4:** This file addresses Criterion 2 in [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561).  
   * **Check** the following in [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561):  
     * Criterion 2  
     * Example 4  

5. `\S0--DESIGN-Inner Trellis Code Stage\RULE_5\Rule_5_LOOP.m`  
   * **Rule 5:** This file addresses the bit assignment described in Remark 5 of [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561).  
   * **Check** the following in [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561):  
     * Appendix C (Table V)  

---

### (III) Design: *Outer-Stage LDPC Code*   [\[S1--DESIGN-Outer LDPC Code Stage\]](https://github.com/arianouri/ISI_Wiretap_Channels.Code_Design/tree/00e9e753840d52aa048bad7ea66d35b805eb785d/%5BSIMULATION_FILES%5D%20Code%20Design/S1--DESIGN-Outer%20LDPC%20Code%20Stage)

1. `\S1--DESIGN-Outer LDPC Code Stage\sub_code_Density_Evolution\main.m`  
   * Given degree distributions, this file evaluates the asymptotic error performance at each iteration using modified density evolution (see Section IV.B of [arXiv:2501.07561 [cs.IT]](https://arxiv.org/abs/2501.07561)).

2. `\S1--DESIGN-Outer LDPC Code Stage\LDPC_ensemble_optimization\main_LP.m`  
   * Given initial degree distributions, this file maximizes the design rate by optimizing the variable-node side degree distributions.

3. 
	* `\S1--DESIGN-Outer LDPC Code Stage\random_code_construction\main_random_PCM.m`
	* `\S1--DESIGN-Outer LDPC Code Stage\random_code_construction\main_encode_PCM.m`
	*  This files are similar to `I.1` and `I.2:`
		`I.1:` `\MC_Encoding_Decoding\main_step1_random_PCM.m`
		`I.2:` `\MC_Encoding_Decoding\main_step2_uppertri_PCM.m`