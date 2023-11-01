#!/usr/bin/env python3

from science.cphmd import theoreticalProtonation

high = 10
low  = 6.5

print('ASP', f"{theoreticalProtonation(high, pKa=3.65):.2f}")
print('ASP', f"{theoreticalProtonation(low,  pKa=3.65):.2f}")

print('GLU', f"{theoreticalProtonation(high, pKa=4.25):.2f}")
print('GLU', f"{theoreticalProtonation(low,  pKa=4.25):.2f}")

print('HIS', f"{theoreticalProtonation(high, pKa=6.53):.2f}")
print('HIS', f"{theoreticalProtonation(low,  pKa=6.53):.2f}")

print('ARG', f"{theoreticalProtonation(high, pKa=13.8):.2f}")
print('ARG', f"{theoreticalProtonation(low,  pKa=13.8):.2f}")

print('LYS', f"{theoreticalProtonation(high, pKa=10.4):.2f}")
print('LYS', f"{theoreticalProtonation(low,  pKa=10.4):.2f}")

print('TYR', f"{theoreticalProtonation(high, pKa=9.84):.2f}")
print('TYR', f"{theoreticalProtonation(low,  pKa=9.84):.2f}")
