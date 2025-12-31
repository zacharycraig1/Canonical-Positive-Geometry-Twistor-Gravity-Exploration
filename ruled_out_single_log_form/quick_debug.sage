load('proper_proof_amplituhedron_hodges.sage')

Z = sample_positive_Z_moment_curve(n=6, seed=0)
twistor = MomentumTwistor(n=6, Z=Z, check_domain=True)

print("Domain OK:", twistor.domain_ok)
print("<0,1>:", twistor.get_angle(0,1))
print("<1,2>:", twistor.get_angle(1,2))
print("<2,0>:", twistor.get_angle(2,0))

H = hodges_6pt_mhv(twistor)
print("Hodges result:", H)

