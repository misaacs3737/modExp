# modExp
Algorithm for fast modular exponentiation modulo exponent-heavy numbers

We present a fast algorithm for modular exponentiation when the factorization of the modulus is known. Let $a,n,m$ be positive integers and suppose $m$ factors canonically as $\prod_{i=1}^k p_i^{e_i}$. Choose integer parameters $t_i\in [0, e_i]$ for $1\le i\le k$. Then we can compute the modular exponentiation $a^n\pmod{m}$ in $O(\max(e_i/t_i)+\sum_{i=1}^k t_i\log p_i)$ steps (i.e. modular operations). We go on to analyze this algorithm mathematically and programmatically, showing significant asymptotic improvement in specific cases.
