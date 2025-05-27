# commented out codeblocks are waiting for new functionalities of tensorflow graph mode
import tensorflow as tf
# from math import pi

from ALP_production.modes.photoproduction import Photoproduction_functions
from ALP_production.general.exotic_constants import alpha_EM#, m_p

from vegasflow import VegasFlow, run_eager, float_me
import logging
logging.getLogger('vegasflow').setLevel(logging.ERROR)
    
def compute_vegas_integral(theta, Ex, x, mx2_4gax2, pt2range, R1, Z_target, photon_density_fun, n_events=1e5, n_iter=10) -> list[float]:
    """returns the photoproduction integral based on the vegas integration method"""
    run_eager() # fixMe: to be removed when tensorflow allows graph evaluation with list comprehension
    tf_setup = Photoproduction_functions_tf(theta, Ex, x, mx2_4gax2, pt2range, [Photoproduction_functions.q2min, (4.49/R1)**2], photon_density_fun, R1, Z_target, Photoproduction_functions.s/Photoproduction_functions.hbarCFm)
    vegas_instance = VegasFlow(2, int(n_events), verbose=False)
    vegas_instance.compile(tf_setup.dsigma_gammaN_intPrep)
    integral, error_est = vegas_instance.run_integration(int(n_iter))
    return integral, error_est

class Photoproduction_functions_tf:
    def __init__(self, theta, Ex, x, mx2_4gax2, pt2range, q2range, photon_denisty_fun, R1, Z, s_hbarCFm):
        import tensorflow as tf
        self.R1 = tf.constant(R1,dtype=tf.float64)
        self.Z = tf.constant(Z,dtype=tf.float64)
        self.theta = tf.constant(theta,dtype=tf.float64)
        self.Ex = tf.constant(Ex,dtype=tf.float64)
        self.mx2_4gax2 = tf.constant(mx2_4gax2,dtype=tf.float64)
        self.s2_2hbarCFm2 = tf.constant(s_hbarCFm**2/2,dtype=tf.float64)
        self.q2min = tf.constant(q2range[0],dtype=tf.float64)
        self.q2max = tf.constant(q2range[1],dtype=tf.float64)
        self.pt2min = tf.constant(pt2range[0],dtype=tf.float64)
        self.dpt2 = tf.constant(pt2range[1] - pt2range[0],dtype=tf.float64)
        self.a_EM = tf.constant(alpha_EM,dtype=tf.float64 )
        # self.pi = tf.constant(pi,dtype=tf.float64)
        # self.mp2 = tf.constant(m_p**2,dtype=tf.float64)
        # self.mup2  = tf.constant(Photoproduction_functions.mu2_p,dtype=tf.float64)
        # self.q2_0 = tf.constant(Photoproduction_functions.q2_0,dtype=tf.float64)
        self.x = x #tf.constant(x,dtype=tf.float64)#
        self.photon_density_fun = photon_denisty_fun
    @staticmethod
    @tf.function
    def custom_3spherical_bessel_j1_x(x):
        """Custom implementation in tf of spherical bessel. 
        approximation error is < 10^8/280."""
        return tf.where(x<float_me(0.01), float_me(1.) - tf.math.pow(x,2)*float_me(0.1), float_me(3)*tf.math.divide(tf.math.divide(tf.math.sin(x),x) - tf.math.cos(x), tf.math.pow(x,2)))
    @tf.function
    def FF(self, q2):
        """Helm Form factor"""
        return  Photoproduction_functions_tf.custom_3spherical_bessel_j1_x(tf.math.sqrt(q2)*self.R1)*tf.exp(-q2 * self.s2_2hbarCFm2)
    @tf.function
    def gamma_density(self, qt2):
        """photon from proton distribution - Budnev """
        qt2_np = qt2.numpy()
        return tf.cast(self.photon_density_fun(self.x, qt2_np), tf.float64) # tf.cast(self.gamma_p_Budnev(qt2), tf.float64)#
    @tf.function
    def dsigma_gammaN_alt(self, pt, phi):
        pE_x = tf.math.pow(pt - self.Ex*self.theta, 2) + float_me(2.)*self.Ex*pt*self.theta*(1. - tf.math.cos(phi))
        return tf.math.divide(self.a_EM * tf.math.pow(self.Ex,2) * pE_x, float_me(4)*tf.math.pow(pE_x+self.mx2_4gax2,2)
                            ) * tf.math.pow(self.Z * self.FF(tf.clip_by_value(pE_x+self.mx2_4gax2, 0, tf.float64.max)), 2)
    @tf.function
    def dsigma_gammaN_intPrep(self, x_array):
        """eq 3.16 in JHEP02(2016)018"""
        phi_scale, pt2_scale = tf.cast(tf.transpose(x_array), dtype=tf.float64)
        pt2 = self.pt2min + pt2_scale*self.dpt2
        pt = tf.math.sqrt(tf.clip_by_value(pt2,0,tf.float64.max))
        phi_lower = tf.math.acos(tf.clip_by_value(( self.mx2_4gax2 + pt2 + tf.math.pow(self.Ex*self.theta,2) - self.q2min)/(float_me(2)*self.Ex*pt*self.theta),-1,1))
        dphi      = tf.math.acos(tf.clip_by_value(( self.mx2_4gax2 + pt2 + tf.math.pow(self.Ex*self.theta,2) - self.q2max)/(float_me(2)*self.Ex*pt*self.theta),-1,1)) - phi_lower
        phi = phi_lower + phi_scale * dphi
        return tf.cast(tf.where(dphi>0., self.dpt2*dphi * self.gamma_density(pt2) *self.dsigma_gammaN_alt(pt, phi), tf.zeros_like(pt)), tf.float64)
    # @tf.function
    # def gamma_p_Budnev(self, qt2):
    #     """photon from proton distribution - Budnev """
    #     Etg2 = qt2 + tf.math.pow(self.x,2)*self.mp2
    #     qi2 = Etg2/(1-self.x)
    #     G2 = 1 / tf.math.pow(1+qi2/Photoproduction_functions.q2_0,4)
    #     C = Photoproduction_functions.mu2_p * G2
    #     D = (float_me(4)*self.mp2*G2 + qi2*Photoproduction_functions.mu2_p*G2)/(float_me(4)*self.mp2+qi2)
    #     return tf.cast(tf.where((Etg2 + self.x) > 1,self.a_EM/self.pi * tf.math.divide(qt2 * (float_me(1)-self.x) * D / Etg2 + tf.math.pow(self.x,2) * C / 2, self.x * Etg2), 0), tf.float64)