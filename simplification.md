That *is* simplify-able, but the main win is to **canonicalize the denominators** first — once you do that, a lot of the repeated “1//24 …” blobs become obviously the same object and you can combine them.

### Step 1: introduce the two detunings

Let
[
\Delta_-=\omega_c-2B,\qquad \Delta_+=\omega_c+2B.
]
Then the four basic reciprocals in your expression become just signs:
[
\frac1{-2B+\omega_c}=\frac1{\Delta_-},\quad
\frac1{2B+\omega_c}=\frac1{\Delta_+},\quad
\frac1{-2B-\omega_c}=-\frac1{\Delta_+},\quad
\frac1{2B-\omega_c}=-\frac1{\Delta_-}.
]

And the “weird quadratic” denominators collapse to perfect squares:
[
-16(B^2)\mp 16B\omega_c-4\omega_c^2=-4(\omega_c\pm 2B)^2=-4\Delta_\pm^2,
]
[
-48(B^2)\mp 48B\omega_c-12\omega_c^2=-12(\omega_c\pm 2B)^2=-12\Delta_\pm^2,
]
[
-48(B^2)+12\omega_c^2 = 12(\omega_c^2-4B^2)=12\Delta_+\Delta_-.
]

### Step 2: name the recurring “difference”

You have this everywhere:
[
\Big(\frac{c_{01}^4 g^4}{-2B+\omega_c}-\frac{c_{01}^4 g^4}{2B+\omega_c}\Big)
= c_{01}^4 g^4\Big(\frac1{\Delta_-}-\frac1{\Delta_+}\Big)
= c_{01}^4 g^4\frac{4B}{\omega_c^2-4B^2}.
]
Call it
[
D := c_{01}^4 g^4\Big(\frac1{\Delta_-}-\frac1{\Delta_+}\Big).
]

### A concrete piece that fully collapses (you can trust this one)

Your quartic “pure creation” term is explicitly
[
\frac{2\Big(\frac{c_{01}^4 g^4}{-2B+\omega_c}-\frac{c_{01}^4 g^4}{2B+\omega_c}\Big)}{-48B^2+12\omega_c^2},a^{\dagger 4}
= \frac{2D}{12(\omega_c^2-4B^2)},a^{\dagger 4}.
]
Since (D=\dfrac{4B c_{01}^4 g^4}{\omega_c^2-4B^2}), this becomes
[
\boxed{;\frac{2B,c_{01}^4 g^4}{3(\omega_c^2-4B^2)^2};a^{\dagger 4};}
]
and the (a^4) term at the end of your expression simplifies to the **same coefficient**:
[
\boxed{;\frac{2B,c_{01}^4 g^4}{3(\omega_c^2-4B^2)^2};a^{4};}
]

So, at minimum, you can replace those two big tails by
[
\frac{2B,c_{01}^4 g^4}{3(\omega_c^2-4B^2)^2},(a^{\dagger 4}+a^4).
]

### Step 3: what the other “square denominators” become

Any occurrence like
[
\frac{D}{-16(B^2)-16B\omega_c-4\omega_c^2}
]
is immediately
[
-\frac{D}{4\Delta_+^2},
]
and similarly
[
\frac{D}{-16(B^2)+16B\omega_c-4\omega_c^2}=-\frac{D}{4\Delta_-^2},
]
[
\frac{D}{-48(B^2)-48B\omega_c-12\omega_c^2}=-\frac{D}{12\Delta_+^2},\quad
\frac{D}{-48(B^2)+48B\omega_c-12\omega_c^2}=-\frac{D}{12\Delta_-^2}.
]

That transformation alone typically shrinks the monster by an order of magnitude, because it turns all those “-16(B^2) ± …” into one of two symbols (\Delta_\pm^2).
