# Units
const °F = u"°F"
const °C = u"°C"
const K = u"K"
const mm = u"mm"
const inch = u"inch"
const cm = u"cm"
const m = u"m"
const s = u"s"
const minute = u"minute"
const hr = u"hr"



# Gelation constants
Ai_yolk, Eai_yolk       = 2.72e50, 3.443e5 # s^-1 and J
Ai_albumen, Eai_albumen = 4.85e60, 4.185e5 # s^-1 and J

R = 8.3144598 # J⋅mol^−1⋅K^−1


# Egg physical parameters as a function of temperature (from literature)
ρ_yolk = T -> 1037.3 - 0.1386 * (T) - 0.0023*(T)^2   # kg/m^3
ρ_albumen = T -> 1043.3 - 0.0115 * (T) - 0.0041*(T)^2  # kg/m^3
c_p_yolk = 3120  # J/(kg·°C)
c_p_albumen = 3800  # J/(kg·°C)
k_yolk = T -> 0.0008*(T) + 0.395 # W/(m·°C)
k_albumen = T -> 0.0013*(T) + 0.5125 # W/(m·°C)