import zlib, base64

"""
This file contains a compiled version of massDynamics.

Do not use this file unless you absolutely cannot get the dynamics to work! If you use this file,
you will have 10 points deducted from your score.

You can use the file by putting the following in your code:

    from massDynamicsCompiled import massDynamics

This `massDynamics` object will then have all of the same methods as in the file `massDynamics.py`:

    * __init__(self, alpha=0.0)
    * update(self, u)
    * f(self, state, F)
    * h(self)
    * rk4_step(self, u)
    * saturate(self, u, limit)

You can then make classes and call the functions as usual. For example:

    from massDynamicsCompiled import massDynamics
    mass = massDynamics()
    y = mass.h()
    
"""


exec(zlib.decompress(base64.b64decode('eJyVUsFuozAQvecrOBpDnZgkq+5K3tPK5xz2VkXISaCxEhNkG6nk69fjAUpbdaVcDPM878145mnT3qxPms60faJc0rQLjZBRzu2UVQbg3eJ4DXEE//SNMvrofi2SU1UnZakb7cuSuOpa5+ranpVYsVUarhOAmPPKV6JpmbJW9eTlJdzu83ju0zHpVfxkz2Pgz5VXYrOlgdTqJX9eUcJZFrUpKRjgVjWnm4kfkj5xlk5SJpTfPkKQRr2J7RhdOLT/kMClAErxCOUQGPwRwl8HNfgY1jd7rMqrNtqLHZNleAKuo2tPYdy4jC4uoRO4BuU7C1dd/llhKmIvm9L5qiUdQD0SzwQCWwV6k/RYpcYCcbW5jGXuIkZM+8qQFVDup5ufgxzANwBnZoCk4AbCl7i8lJKnYQ/0PvwV9E7pOsMEin6BaTndkHfDpJnE/AMF0XSf72eNQ13s/Rw5Y8+TQ2eN9x8a/KgzDGCa1GzQkqMcDgdV8zhJWXy9yYa1Lgsq+ZC2/n9aMaRtvk+jco1JsysxSvygRPIMhOBYZ3KT4msmb+BrcnQFvEnXiTq44IffEQMoGCr+4wZeGzTLMJ5u8Q/Gq1YN')))
