# ThreeBody
A n-body problem solver with order-4 Runge-Kutta (RK4) method.

## Math Model
Assume there are $n$ stars in $\mathbb{R}^3$ space, with mass $m_i$, position $s_i$ and and velocity $v_i$ respectively.

By Newton's laws, we have
$$\frac{\mathrm{d}s_i}{\mathrm{d}t}=v_i,$$
$$\frac{\mathrm{d}v_i}{\mathrm{d}t}=\mathrm{G}\cdot\sum_{j\not=i}\frac{s_j-s_i}{\lvert s_j-s_i \rvert}\cdot m_j$$

Discretize the differential equations
$$s^{n+1}_i(v^{(n)}_i) = s_i^{(n)} + v^{(n)}_i \cdot \Delta t$$
$$v^{(n+1)}_i(s^{(n)}) = v_i^{(n)} +  \mathrm{G}\cdot\sum_{j\not=i}\frac{s^{(n)}_j-s^{(n)}_i}{\lvert s^{(n)}_j-s^{(n)}_i \rvert}\cdot m_j\cdot \Delta t.$$

Let
$$f(v) = v,$$
$$g(s, i) = \mathrm{G}\cdot\sum_{j\not=i}\frac{s_j-s_i}{\lvert s_j-s_i \rvert}\cdot m_j.$$
By order-4 Runge-Kutta method, we have
$$k^{(v_i)}_1 = f(v_i^{(n)}),$$
$$k^{(s_i)}_1 = g(s_i^{(n)}, i);$$
$$k^{(v_i)}_2 = f(v_i^{(n)}+k^{(v_i)}_1 \cdot\frac{\Delta t}{2}),$$
$$k^{(s_i)}_2 = g(s_i^{(n)} + k^{(v_i)}_1 \cdot \frac{\Delta t}{2}, i);$$
$$k^{(v_i)}_3 = f(v_i^{(n)}+k^{(v_i)}_2 \cdot\frac{\Delta t}{2}),$$
$$k^{(s_i)}_3 = g(s_i^{(n)} + k^{(v_i)}_2 \cdot \frac{\Delta t}{2}, i);$$
$$k^{(v_i)}_4 = f(v_i^{(n)}+k^{(v_i)}_3 \cdot\Delta t),$$
$$k^{(s_i)}_4 = g(s_i^{(n)} + k^{(v_i)}_3 \cdot \Delta t, i).$$
Finally we get the approximate position and velocity at $t + \Delta t$
$$v^{(n+1)}_i = \frac{\Delta t}{6}(k^{(v_i)}_1 + 2k^{(v_i)}_2 + k^{(v_i)}_3 + k^{(v_i)}_4)$$
$$s^{(n+1)}_i = \frac{\Delta t}{6}(k^{(s_i)}_1 + 2k^{(s_i)}_2 + k^{(s_i)}_3 + k^{(s_i)}_4)$$


## Related Material
OJ problem and solution motivation: [洛谷 P3945 三体问题](https://www.luogu.com.cn/problem/P3945)

ODE numerical solution methods: [Runge-Kutta Methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
