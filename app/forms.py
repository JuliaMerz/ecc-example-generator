from wtforms import Form, BooleanField, IntegerField, validators

class ECCForm(Form):
    n = IntegerField('Length of original message (n in notes)', [validators.Required()])
    k = IntegerField('Number of errors to protect against (k in notes):', [validators.Required()])
    mod = IntegerField('Prime number we are working modulo (must be prime and bigger than n + 2k)\nEnter \"0\" to default to next highest prime:', [validators.Required()])
    show_inv = BooleanField('Multiplicative inverse / EGCD?')
    show_rref = BooleanField('Row Reduction to solve system of equations?')
    show_lagrange = BooleanField('Lagrange Polynomial Interpolation?')
    show_orig = BooleanField('The original polynomial/message at the beginning? Leave unchecked for practice problems.')
    show_bw = BooleanField('The Berlekamp-Welch algorithm to construct a system of equations?')
