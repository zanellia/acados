#ifndef {{ ocp.con_p.name }}_P_CONSTRAINT
#define {{ ocp.con_p.name }}_P_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if ocp.dims.npd > 0 %}
// implicit ODE
int {{ ocp.con_p.name }}_p_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.con_p.name }}_p_constraint_work(int *, int *, int *, int *);
const int *{{ ocp.con_p.name }}_p_constraint_sparsity_in(int);
const int *{{ ocp.con_p.name }}_p_constraint_sparsity_out(int);
int {{ ocp.con_p.name }}_p_constraint_n_in();
int {{ ocp.con_p.name }}_p_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ocp.con_p.name }}_P_CONSTRAINT
