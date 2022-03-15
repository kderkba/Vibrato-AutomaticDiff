# Vibrato and AD for Option Greeks

(Preview)

Vibrato method and the Vibrato + Automatic Differentiation from the paper 'Vibrato and automatic differentiation for high order derivatives and sensitivities' https://arxiv.org/abs/1606.06143 are implemented

The following results are obtained from the algorithm:

## European Options
![Gamma European Option](https://github.com/kderkba/Vibrato-AutomaticDiff/blob/main/VAD_gamma_euro.png)

![Convergence](https://github.com/kderkba/Vibrato-AutomaticDiff/blob/main/convergence_gamma_VAD.png)


## Digital Options
![Vega Digital Option](https://github.com/kderkba/Vibrato-AutomaticDiff/blob/main/VAD_vega_digital.png)

As shown in the paper, the algorithm can be applied to different types of options. 
Futher development will include the implementation of Basket option greeks and Heston greeks.
