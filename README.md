# Black-to-White hole transition

## Structure of computation

The code computes the Black-to-White hole transitin amplitude:

$$
W_{\alpha} (j_{\pm}, j_{0}, \zeta_{\pm}, \zeta_0; \Delta l) = \sum_{ j_{ab}^{\pm}, j_{a} } w_{\alpha} \left( \sum_{i_5} d_{i_5} \prod_{\pm} V_{\gamma,  n_{f}^{\pm}}^{coh} \left(j_{ab}^{\pm}, j_{a}, i_5; \Delta l \right) \right) 
$$

and the crossing time:

$$
\tau_{\alpha} \left( m \right) = \int_{0}^{\infty} dT \  T \ P_{\alpha} \left( m | T \right)  \ ,
$$

where:

$$
P_{\alpha} \left( m | T \right) = \frac{\mu_{\alpha} (m, T) \ |W_{\alpha} (m , T )|^2}{\int_{0}^{\infty} dT \ \mu_{\alpha} (m, T) \ |W_{\alpha} (m, T)|^2} \ .
$$

The computation is divided into multiple steps:

- *bin/vertex_computation.jl* computes all the required full vertex amplitude tensors. A supercluster is recommended for high values of $\Delta l$.

- *bin/vertex_contraction.jl* evalueates the term: $\left( \sum_{i_5} d_{i_5} \prod_{\pm} V_{\gamma,  n_{f}^{\pm}}^{coh} \left(j_{ab}^{\pm}, j_{a}, i_5; \Delta l \right) \right)$ for each value of $j_{ab}^{\pm}, j_{a}$.

- *bin/weight_factor.jl* computes the weight factor $w_{\alpha}$ (whose expression is not explicitly reported here) for each value of $j_{ab}^{\pm}, j_{a}$ and $T_1 \dots T_N$. It requires the number $N$ of T samples as additional parameter.

- *bin/amplitude_assembling.jl* puts together all the terms and computes the total amplitude.

The data are plotted in the *analysis.ipynb* notebook, inside the folder `notebooks/Julia`.

#TODO ...

**The Julia codes are parallelized on the available cores.** It is therefore advisable for the performance to parallelize the codes keeping into account the number of CPU cores present on the system.

A full list of the employed Julia packages can be found in `./inc/pkgs.jl`. **Before executing the source codes, all packages must be installed.**

**The Julia's Just-in-Time compiler is such that the first execution of functions is considerably slower that following ones, and it also allocates much more memory**. To avoid this, you can use the [DaemonMode package](https://github.com/dmolina/DaemonMode.jl).

## Usage

To run the Julia scripts (on a single machine with the synthax below) you can use the following command:

```
$JULIA_EXECUTABLE_PATH   -p   [N-1]    $JULIA_CODE_PATH   $STORE_FOLDER    $DATA_SL2CFOAM_FOLDER
```

where [N-1] is the number of workers and:

- `STORE_FOLDER`: folder where data are saved

- `DATA_SL2CFOAM_FOLDER`: folder with fastwigxj tables where booster functions are stored/retrieved.

