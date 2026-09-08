/***********************************************/
/**
* @file kernelGeoid.h
*
* @brief Integral kernel of geoid computation.
* (= Poisson Kern * gamma).
* @see Kernel
*
* @author Torsten Mayer-Guerr
* @date 2003-09-20
*
*/
/***********************************************/

#ifndef __GROOPS_KERNELGEOID__
#define __GROOPS_KERNELGEOID__

// Latex documentation
#ifdef DOCSTRING_Kernel
static const char *docstringKernelGeoid = R"(
\subsection{GeoidHeight}\label{kernelType:geoidHeight}
The geoid height is defined by Bruns formula
\begin{equation}
N = \frac{1}{\gamma}T
\end{equation}
with $T$ the disturbance potential and the normal gravity
\begin{equation}\label{normalgravity}
  \gamma  = \gamma_0 - 0.30877\cdot 10^{-5}/s^2(1-0.00142\sin^2(B))h + 0.72\cdot 10^{-12}/(m\,s^2)h^2
\end{equation}
and
\begin{equation}
  \gamma_0 = \frac{a\gamma_a\cos^2(B)+b\gamma_b\sin^2(B)}{\sqrt{a^2\cos^2(B)+b^2\sin^2(B)}}
\end{equation}
where $h$ is the ellipsoidal height in meter, $B$ the ellipsoidal latitude,
$a$, $b$ the semi-axes of the GRS80 ellipsoid and $\gamma_a=9.7803267715\,m/s^2$,
$\gamma_b=9.8321863685\,m/s^2$ the normal gravity at the equator and at the
pole, as given in GRS80 (Moritz 1980).

The kernel is given by
\begin{equation}
K(\cos\psi,r,R) = \gamma\frac{R(r^2-R^2)}{l^3},
\end{equation}
and the coefficients in \eqref{eq.kernel} are
\begin{equation}
k_n = \gamma.
\end{equation}
)";
#endif

/***********************************************/

#include "classes/kernel/kernel.h"

/***** CLASS ***********************************/

/** @brief Integral kernel of geoid computation.
* @ingroup kernelGroup
* (= Poisson Kern * gamma).
* @see Kernel */
class KernelGeoid : public Kernel
{
public:
  KernelGeoid(Config &/*config*/) {}
  Double   kernel             (const Vector3d &p, const Vector3d &q) const;
  Double   radialDerivative   (const Vector3d &p, const Vector3d &q) const;
  Vector3d gradient           (const Vector3d &p, const Vector3d &q) const;
  Tensor3d gradientGradient   (const Vector3d &p, const Vector3d &q) const;
  Double   inverseKernel      (const Vector3d &p, const Vector3d &q, const Kernel &kernel) const;
  Double   inverseKernel      (const Time &time, const Vector3d &p, const GravityfieldBase &field) const;
  Vector   coefficients       (const Vector3d &p, UInt degree) const;
  Vector   inverseCoefficients(const Vector3d &p, UInt degree, Bool interior) const;
};

/***********************************************/

#endif
