#ifndef LIMITERS_HPP
#define LIMITERS_HPP

__device__ void connectivity_stats(double x_i, double y_i, double nx, double ny, double power, double conn_x, double conn_y, double& sig_del_x_sqr, double& sig_del_y_sqr, double& sig_del_x_del_y, double& delta_x, double& delta_y, double& delta_s_weights, double& delta_n_weights);

__device__ void calculate_qtile(double qtilde_i[4], double qtilde_k[4], Point* globaldata, int idx, int conn, double delta_x, double delta_y, double vl_const, double gamma, int limiter_flag, double phi_i[4], double phi_k[4]);

__device__ void venkat_limiter(double qtilde[4], double vl_const, Point* globaldata, int index, double gamma, double phi[4]);

__device__ void VLBroadcaster(double q[4], double qtilde[4], double max_q[4], double min_q[4], double phi[4], double epsi, double del_pos, double del_neg);

__device__ void qtilde_to_primitive(double result[4], double qtilde[4], double gamma);

__device__ void update_delf(double sig_del_x_del_f[4], double sig_del_y_del_f[4], double G_k[4], double G_i[4], double delta_s_weights, double delta_n_weights);

template <class Type>
bool isNan(Type var);

#endif
