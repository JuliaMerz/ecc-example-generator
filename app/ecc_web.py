import random
import operator
import collections
import itertools
from itertools import izip
from math import fabs


"""
Practice problem generator with solutions for Reed Solomon codes as presented in CS70, Spring 2013 at UC Berkeley.
Written by Quinn Johnson.

Made web compatible Spring 2014. If you wish to use this code for anything else, it is 
recommended you find the original, since changes made to allow use for websites disallow 
for independent use of any function besides run_web.
Modified by Sebastian Merz
"""


#String wrapper class to allow us to access a string within functions.
class CludgeString:
    def __init__(self):
        self.text = ''
    def clear(self):
        self.text = ''
    def append(self, string):
        self.text = self.text + str(string)
    def append_line(self, string):
        self.text = self.text + str(string) + "\n"
    def print_output(self):
        return self.text

output = CludgeString()
#EGCD as presented in the notes. Has extra parameter of "show" to show work.
def egcd(x, y,show):
    if y == 0:
        if show:
            output.append_line("x = " + str(x) + ", y = " + str(y))
            output.append_line("d = " + str(x) + ", a = 1, b = 0")
        return (x, 1, 0)
    else:
        if show:
            output.append_line("x = " + str(x) + ", y = " + str(y))
        (d, a, b) = egcd(y, x % y,show)
        if show:
            temp = (a - (x / y) * b)
            output.append_line("d = " + str(d) + ", a = " + str(b) + ", b = " + str(temp))
    return ((d , b , a - (x / y) * b))


#Function that calculates the multiplicative inverse of a number by calling egcd. 
#Has parameter of show to show work.
def multinv(num, mod,show):
    if num < 0:
        num = num%mod

    x = egcd(mod, num,show)[2]
    if x < 0:
        if show:
            output.append_line(str(x) + " = " + str(mod + x) + " mod " + str(mod))
        x = mod + x
    return x


#Function that takes in a matrix and returns a row-reduced version, modulo a number mod
#Has show parameter to show work
def rowreduce_mod(matrix,mod,show):
    rows = len(matrix)
    cols = len(matrix[0])
    i,j = 0,0
    
    #Loop until reduced
    while True:
        #If out of bounds, reduced
        if i >= rows or j >= cols:
            break
        
        #Find a matrix that isn't zero (contains pivot)
        not_zero_row = i
        if matrix[i][j] == 0:
            while not_zero_row < rows and matrix[not_zero_row][j] == 0:
                not_zero_row += 1
        
        #Jumps to next column if no pivot is found
        if not_zero_row == rows:
            j += 1
            continue
        
        #Swaps pivot row up
        temp = matrix[i]
        matrix[i] = matrix[not_zero_row]
        matrix[not_zero_row] = temp

        #Normalizes pivot row so that pivot = 1
        pivot = matrix[i][j]
        piv_inv = multinv(pivot,mod,False)
        matrix[i] = [((x * piv_inv) % mod) for x in matrix[i]]

        #Adds multiple of pivot row to create 0's above and below
        for row in range(0, rows):
            if row == i:
                continue
            if matrix[row][j] != 0:
                matrix[row] = [((y - matrix[row][j]*x)%mod) for (x,y) in zip(matrix[i], matrix[row])]
        
        i += 1
        j += 1
        
        #Option added to show intermediary steps
        if show:
            output.append_line("")
            for row in matrix:
                output.append_line(row)
            output.append_line("")
            pause()
    return matrix


#Taken from http://stackoverflow.com/questions/5413158/multiplying-polynomials-in-python
#Fantastic implementation directly in python. 
#I opted to use it instead of a library to make the entire operations in this program transparent
#Methods added for our purposes, including modular arithmetic and division
#Some of the original methods were also modified to fit the needs of this program

class Polynomial(object):
    def __init__(self, *args):
        """
        Create a polynomial in one of three ways:

        p = Polynomial(poly)           # copy constructor
        p = Polynomial([1,2,3 ...])    # from sequence
        p = Polynomial(1, 2, 3 ...)    # from scalars
        """
        super(Polynomial,self).__init__()
        if len(args)==1:
            val = args[0]
            if isinstance(val, Polynomial):                # copy constructor
                self.coeffs = val.coeffs[:]
            elif isinstance(val, collections.Iterable):    # from sequence
                self.coeffs = list(val)
            else:                                          # from single scalar
                self.coeffs = [val+0]
        else:                                              # multiple scalars
            self.coeffs = [i+0 for i in args]
        #self.trim()

    def __add__(self, val):
        "Return self+val"
        if isinstance(val, Polynomial):                    # add Polynomial
            res = [a+b for a,b in itertools.izip_longest(self.coeffs, val.coeffs, fillvalue=0)]
        else:                                              # add scalar
            if self.coeffs:
                res = self.coeffs[:]
                res[0] += val
            else:
                res = val
        return self.__class__(res)

    def __call__(self, val):
        "Evaluate at X==val"
        res = 0
        pwr = 1
        for co in self.coeffs:
            res += co*pwr
            pwr *= val
        return res

    def __eq__(self, val):
        "Test self==val"
        if isinstance(val, Polynomial):
            return self.coeffs == val.coeffs
        else:
            return len(self.coeffs)==1 and self.coeffs[0]==val

    def __mul__(self, val):
        "Return self*val"
        if isinstance(val, Polynomial):
            _s = self.coeffs
            _v = val.coeffs
            res = [0]*(len(_s)+len(_v)-1)
            for selfpow,selfco in enumerate(_s):
                for valpow,valco in enumerate(_v):
                    res[selfpow+valpow] += selfco*valco
        else:
            res = [co*val for co in self.coeffs]
        return self.__class__(res)

    def __neg__(self):
        "Return -self"
        return self.__class__([-co for co in self.coeffs])

    def __pow__(self, y, z=None):
        raise NotImplemented()

    def _radd__(self, val):
        "Return val+self"
        return self+val

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, self.coeffs)

    def __rmul__(self, val):
        "Return val*self"
        return self*val

    def __rsub__(self, val):
        "Return val-self"
        return -self + val

    def __str__(self):
        "Return string formatted as aX^3 + bX^2 + c^X + d"
        res = []
        for po,co in enumerate(self.coeffs):
            if co:
                if po==0:
                    po = ''
                elif po==1:
                    po = 'X'
                else:
                    po = 'X^'+str(po)
                res.append(str(co)+po)
        if res:
            res.reverse()
            return ' + '.join(res)
        else:
            return "0"

    def __sub__(self, val):
        "Return self-val"
        return self.__add__(-val)

    def trim(self):
        "Remove trailing 0-coefficients"
        _co = self.coeffs
        if _co:
            offs = len(_co)-1
            if _co[offs]==0:
                offs -= 1
                while offs >= 0 and _co[offs]==0:
                    offs -= 1
                del _co[offs+1:]
    def mod(self, num):
        for i in range(len(self.coeffs)):
            self.coeffs[i] = int(self.coeffs[i])%num
        return self

    def degree(self, poly):
        while poly and poly[-1] == 0:
            poly.pop()   # normalize
        return len(poly)-1
 
    def poly_div(self, N, D,mod):
        dD = self.degree(D)
        dN = self.degree(N)
        if dD < 0: raise ZeroDivisionError
        if dN >= dD:
            q = [0] * dN
            while dN >= dD:
                d = [0]*(dN - dD) + D
                mult = q[dN - dD] = N[-1] * multinv((d[-1]),mod,False)
                d = [coeff*mult for coeff in d]
                N = [fabs ( (coeffN - coeffd)%mod ) for coeffN, coeffd in izip(N, d)]
                dN = self.degree(N)
            r = N
        else:
            q = [0]
            r = N
        return q, r

    def div(self,x,mod):
        a = self.coeffs
        b = x.coeffs
        c = self.poly_div(a,b,mod)
        d = [int(x) for x in c[0]]
        return Polynomial(d)


#Lagrange Interpolation
#You pass in a list of points in the form [[x1,y1],[x2,y2],...], a mod, and a show parameter
#Returns a Polynomial that matches all points, of degree 1 less than number of points.
#Show parameter is there to allow showing of work

def lagrange(points,mod,show):
    #Construct list of Delta polynomials
    delta_list = []
    #Optional statement for showing work
    if show:
        output.append_line("Sample points are: " + str(points))

    #Go through each point, creating a delta based on the x value
    #Strings are constructed during the process to allow shown work
    for i in range(len(points)):
        numer = Polynomial(1)
        denom = 1
        numer_string = ""
        denom_string = ""
        for j in range(len(points)):
            if i != j:
                numer_string = numer_string + "(" + str(Polynomial((-1*points[j][0]),1)) + ") "
                denom_string = denom_string + "(" + str(points[i][0]) + " - " + str(points[j][0]) + ") "
                numer = numer * Polynomial((-1 * points[j][0]),1)
                denom = denom * (points[i][0] - points[j][0])
        denom = multinv(denom, mod,False)
        delta_poly = numer * denom
        delta_list.append(delta_poly)
        #Optional statement for showing work
        if show:
            output.append_line("Delta " + str(points[i][0]) + ": " + numer_string + " * inverse of " + denom_string + " = " + str(delta_poly) + " = " + str(delta_poly.mod(mod)) + " mod " + str(mod))
    val = Polynomial(0)
    if show:
        output.append_line("\nWe now add these all together, weighting them with the appropriate Y value:")
    final_string = ""

    #Combining weighted deltas to get the final polynomial
    for i in range(len(points)):
        final_string = final_string + "(" + str(points[i][1]) + ")(" + str(delta_list[i]) + ")"
        if i < len(points) - 1:
            final_string = final_string + " + "
        val = val + (delta_list[i] * points[i][1])
    if show:
        output.append_line(final_string)
    return val.mod(mod)


#Solver that takes in a matrix, modulus, and show parameter and solves it.
#For simplicity, makes the assumption that any free variables equal 1
#Returns a vector of coeffecients
def matrix_solver(mat,mod,show):
    matrix = rowreduce_mod(mat,mod,show)
    sol_vec = []
    i,j = 0,0
    rows = len(matrix)
    cols = len(matrix[0])-1

    while True:
        if i >= rows or j >= cols:
            break
        
        not_zero_row = i
        if matrix[i][j] != 1:
            while not_zero_row < rows and matrix[not_zero_row][j] == 0:
                not_zero_row += 1
        
        if not_zero_row == rows:
            j += 1
            sol_vec.append(1)
            continue
        
        temp = matrix[not_zero_row][cols]
        for inc in range(cols):
            if inc != j:
                temp = temp - matrix[not_zero_row][inc]
        temp = temp%mod
        sol_vec.append(temp)
        i += 1
        j += 1

    return sol_vec


#Simple check to see if a number n is prime
def prime(n):
    if n % 2 == 0 or n % 3 == 0:
        return False
    for f in range(5, int(n ** .5),6):
        if n % f == 0 or n % (f + 2) == 0:
            return False
    return True


#Function to construct system of equations using the Berlekamp-Welch Decoding Algorithm
def bw_alg(length, rec, q_deg, e_deg, mod, show):
    p_list = []
    #output.append_line(the form of Q(x) and E(x) if user wants)
    if show:
        output.append_line("Because Q(x) is of degree " + str(q_deg) + ", and E(x) is of degree " + str(e_deg) + ", he knows:")
        q_temp = q_deg
        q_string = ""
        while not q_temp == 0:
            q_string += "(a" + str(q_temp) + ")X^(" + str(q_temp) + ") + "
            q_temp = q_temp - 1
        q_string += "(a0) mod " + str(mod)
        
        e_temp = e_deg
        e_string = "X^(" + str(e_temp) + ") + "
        e_temp = e_temp - 1
        while not e_temp == 0:
            e_string += "(b" + str(e_temp) + ")X^(" + str(e_temp) + ") + "
            e_temp = e_temp - 1
        e_string += "(b0) mod " + str(mod)
        output.append_line("R(x) is defined by the points Bob received.")
        output.append_line("Q(x) = " + q_string)
        output.append_line("E(x) = " + e_string + "\n")

    #For each value of x, calculate coefficients
    for i in range(1,length+1):
        qvec = []
        evec = []
        for j in range(q_deg+1):
            qvec.append(pow(i,j)%mod)
        for k in range(e_deg+1):
            evec.append(pow(i,k)%mod)
        
        if show:
            output.append_line("Equation " + str(i) + " for x = " + str(i) + ":")
            num = "(" + str(i) + ")"
            output.append_line("Q" + num + "\t= " + q_string.replace('X', num))
            q_poly = str(Polynomial(qvec)).replace('X^', 'a') + "a0"
            q_poly = q_poly.replace('X', 'a1')
            output.append_line("\t= " + q_poly + " mod " + str(mod))
            output.append_line("E" + num + "\t= " + e_string.replace('X', num))
            e_poly = str(Polynomial(evec)).replace('X^' + str(e_deg), "")
            e_poly = e_poly.replace('X^', 'b') + "b0"
            e_poly = e_poly.replace('X', 'b1')
            output.append_line("\t= " + e_poly + " mod " + str(mod))
            output.append_line("Q" + num + "=R" + num + "E" + num + ":\t" + q_poly + " = " + str(rec[i-1]) + "(" + e_poly + ") mod " + str(mod))

        for l in range(len(evec)):
            evec[l] = (evec[l] * rec[i-1])%mod
        
        if show:
            e_poly = str(Polynomial(evec)).replace('X^' + str(e_deg), "")
            e_poly = e_poly.replace('X^', 'b') + "b0"
            e_poly = e_poly.replace('X', 'b1')
            output.append_line("\t=>\t" + q_poly + " = " + e_poly + " mod " + str(mod))

        qvec.reverse()
        evec.reverse()
        #Combine into one vector
        vec = qvec
        if show:
            q_temp = q_deg
            final_eq = ""
            for x in vec:
                final_eq += str(x) + "a" + str(q_temp) + " + "
                q_temp = q_temp - 1
            e_temp = e_deg - 1
        for j in range(1,e_deg+1):
            app_val = ((-1)*evec[j])%mod
            vec.append(app_val)
            if show:
                final_eq += str(app_val) + "b" + str(e_temp)
                if not e_temp == 0:
                    final_eq += " + "
                e_temp = e_temp - 1

        vec.append(evec[0])
        
        if show:
            final_eq += " = " + str(evec[0]) + " mod " + str(mod)
            output.append_line("\t=>\t" + final_eq + "\n\n")
            pause()
        p_list.append(vec)
    return p_list


#Gets intial values from the user
def get_user_input():
    try:
        output.append_line("\nPlease enter variables relevant to the Reed-Solomon problem:")

        n = int(raw_input("\nLength of original message (n in notes): \n"))
        k = int(raw_input("\nNumber of errors to protect against (k in notes): \n"))
        mod = int(raw_input("\nPrime number we are working modulo (must be prime and bigger than n + 2k)\nEnter \"0\" to default to next highest prime: \n"))
        if mod == 0:
            i = n + 2*k + 1
            while not prime(i):
                i += 1
            mod = i
            output.append_line("Calculations will be done modulo " + str(mod) + "\n")
        
        output.append_line("\nNow enter what steps you want to see covered in detail. All set to zero provides a minimal solution.")
        show_inv = str(raw_input("\nMultiplicative inverse / EGCD? (y/n): \n"))
        if show_inv == 'n' or show_inv == 'N':
            show_inv = False
        
        show_rref = str(raw_input("\nRow Reduction to solve system of equations? (y/n): \n"))
        if show_rref == 'n' or show_rref == 'N':
            show_rref = False
        
        show_lagrange = str(raw_input("\nLagrange Polynomial Interpolation? (y/n): \n"))
        if show_lagrange == 'n' or show_lagrange == 'N':
            show_lagrange = False

        show_orig = str(raw_input("\nThe original polynomial/message at the beginning? Choose no for practice problems. (y/n): \n"))
        if show_orig == 'n' or show_orig == 'N':
            show_orig = False

        show_bw = str(raw_input("\nThe Berlekamp-Welch algorithm to construct a system of equations? (y/n): \n"))
        if show_bw == 'n' or show_bw == 'N':
            show_bw = False

        output.append_line("\n\n"    )
        
    except:
        output.append_line("Input error")
        exit()
    return [n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw]


#Pauses for user input to continue
def pause():
    #Stop pausing from slowing us down
    #raw_input("Press Enter to continue...")
    output.append_line("")


#Main interactive function
def run_web(n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw):
    print n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw
    random.seed()
    if mod == 0:
        i = n + 2*k + 1
        while not prime(i):
            i += 1
        mod = i
        output.append_line("Calculations will be done modulo " + str(mod) + "\n")
        

    #We use this to make conversion from js possible without changing architecture completely
    #We then replace all print with output.append_line
    output.clear()
    
    #Get values from user
    #We remove this so we can get input from web
    #[n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw] = get_user_input()
    length = n + (k*2)
    orig = [0]*n
    rec = [0]*length
    point_list = []
    
    #Generate original message
    for i in range(n):
        orig[i] = random.randint(0,mod-1)
        point_list.append([i+1,orig[i]])

    output.append_line("For this problem, Alice wants to send a message of length " + str(n) + " to Bob.")
    output.append_line("She also wants to ensure that it can survive up to " + str(k) + " general errors.")
    output.append_line("She realizes that through the Reed-Solomon scheme, she can do so with a message of length " + str(length) + " over GF(" + str(mod) +").")
    if (show_orig):
        output.append_line("\nHer message is:\n " + str(orig) + "\n")
    
    #Calculating degrees of polynomials
    q_deg = n + k -1
    e_deg = k
    p_deg = n - 1

    output.append_line("We know that if we create an original polynomial P(x) of degree " + str(p_deg) + " passing through the points of her message, we can extend it.")
    output.append_line("From the extended (and possibly corrupted) message, we can obtain P(x) from two polynomials:")
    output.append_line("The error polynomial, E(x), of degree " + str(e_deg) + " and P(x)E(x) = Q(x), of degree " + str(q_deg) + ".")

    #Many optional output.append_line(statements if user wanted to see original information before being sent)
    if (show_orig and show_lagrange):
        output.append_line("\nHere is the interpolation she does to find her original polynomial, P(x)\n")
    
    show_lagrange_early = show_lagrange and show_orig
    #Calculate the original polynomial
    orig_poly = lagrange(point_list,mod,show_lagrange_early)
    if (show_orig):
        output.append_line("\nHer original polynomial, obtained through interpolation, is: \n" + str(orig_poly))
    #Creates the message Alice sends
    for i in range(length):
        rec[i] = orig_poly(i+1)%mod

    if (show_orig):
        output.append_line("Using this, she sends out the full message: " + str(rec) + "\n")
    #Randomizes the message to create errors
    for _ in range(k):
        a = random.randint(0,length-1)
        rec[a] = random.randint(0,mod-1)

    pause()
    #Set of output.append_line(statements and computations if user wants to see Multiplicative Inverse calculations)
    if show_inv:
        output.append_line("Because you, the user, asked to see the computation of the inverses, we will do so for all values less than our modulus.\n")
        inv_list = [-1]*mod
        for i in range(1,mod):
            if inv_list[i] == -1:
                output.append_line("We have not yet calculated the inverse of " + str(i) + " so we will do that with egcd:\n")
                a = multinv(i, mod, show_inv)
                inv_list[a] = i
                inv_list[i] = a
                output.append_line("\nFrom this, we know that " + str(i) + " is " + str(a) + "\'s inverse, and vice versa.\n")
                pause()
            else:
                output.append_line("We have already calculated the inverse of " + str(i) + " so we will skip it. \n")
                pause()

    #Explain how to construct system of equations:
    #Could possibly be done with the """ """ construct, but that may interfere with doc-strings if I add them.
    output.append_line("Bob receives the message: \n " + str(rec) + "\nwhich he knows may be wrong in up to " + str(k) + " places.\n")
    output.append_line("With this, he seeks to create a system of equations of the form Q(x) = R(x)E(x) for each point he received.")
    output.append_line("Using this, he can solve to reconstruct the original P(x).")

    if show_bw:
        output.append_line("He knows that Q(x) has " + str(q_deg + 1) + " degrees of freedom, and E(x) has " + str(e_deg) + " degrees of freedom (because the leading coeffecient is always 1).\n")
        output.append_line("To do so, he evaluates Q(x) = R(x)E(x) for each value of x he has received.")
        output.append_line("From his received message, R(x) corresponds to a received value. For instance, R(1) = " + str(rec[0]) + ".")
        output.append_line("We can construct Q(x) and E(x) as follows:\n")

        output.append_line("Q(x) = a_(n+k-1) * x^(n+k-1) + ... + a_1 * x + a_0")
        output.append_line("E(x) = x^k + b_(k-1) * x^(k-1) + ... + b_1 * x + b_0")
        output.append_line("Note here that the leading coeffecient is always 1.")
        output.append_line("Evaluating at each point, he constructs the following equations:\n\n")
        pause()    
    poly_list = bw_alg(length, rec, q_deg, e_deg, mod, show_bw)
    
    #Display system of equations
    output.append_line("By plugging in values for each x he was given, he creates a system of equations.")
    output.append_line("Because it is a system of " + str(q_deg + e_deg - 1) + " linear equations in " + str(q_deg + e_deg - 1) + "variables, he represents it as the following matrix:")
    output.append_line("(Note: If you have not taken Math 54 or are unfamiliar with the row-reduction method of solving systems of equations, I advise you to look it up online. It will make your life much easier.)\n\n")
    for row in poly_list:
        output.append_line(row)
    
    pause()
    output.append_line("\nHe then row-reduces this matrix to solve for each unknown: ")

    #Solution taken is one where all free variables are 1
    sols = matrix_solver(poly_list,mod,show_rref)
    output.append_line("\nAfter row reduction, he arrives at the following set of values (note that he takes all free variables to be 1): ")
    output.append_line(str(sols) + "\n")

    #Constructing and displaying solutions as 2 sets of coefficients
    qvec = sols[:q_deg+1]
    output.append_line("The first portion of those values correspond to the coeffecients of Q(x), from highest degree to lowest: \n" +  str(qvec) + "\n")
    qvec.reverse()
    evec = sols[q_deg+1:]
    output.append_line("The later portion of those values correspond to the coeffecients of E(x), from highest degree to lowest: \n" + str(evec) + "\n")
    output.append_line("Note that for the coeffecients of E(x), it is understood that the leading coefficient is 1, and is not in the above list.\n")
    evec.reverse()
    evec.append(1)
    
    #Creating polynomials
    qx = Polynomial(qvec)
    ex = Polynomial(evec)
    
    pause()
    output.append_line("From this, he can reconstruct the Q(x) and E(x): ")
    output.append_line("Q(x) = " + str(qx))
    output.append_line("E(x) = " + str(ex) + "\n")
    
    #Calculating original polynomial
    pause()
    px = (qx.div(ex,mod)).mod(mod)
    output.append_line("Bob then divides Q(x) by E(x) to finally get back to the original polynomial P(x):\n" + str(px))
    output.append_line("\nBecause the original message lies along the values of P(x), Bob now has sufficient information to find any errors in his received copy.")

    #Finding roots of E(x)
    output.append_line("Bob also knows that the roots of E(X) correspond to the locations of errors, so he solves:")
    output.append_line(str(ex) + " = 0")
    output.append_line("and finds that the message was corrupted at positions:\n")
        
    roots = []
    for i in range(0,mod):
        val = ex(i)%mod
        if val == 0:
            roots.append(i)
    
    #If no errors, special message
    if len(roots) == 0:
        output.append_line("Huh. There were no errors. I guess Bob just got lucky today, didn't he?")
    #Otherwise, output.append_line(out original)
    else:
        pause()
        err_locs = ""
        for err in roots:
            err_locs += str(err) + "\t"
        output.append_line(err_locs + "\n")
        output.append_line("So, he recalculates to find that the actual values are:\n")
        for err in roots:
            output.append_line("Position " + str(err) + ": " + str(px(err)%mod) + ", not " + str(rec[err-1]))
    
        output.append_line("\nFinally, he fixes the errors he found to arrive back at the original transmitted message:")
        fixed = []
        for i in range(length):
            fixed.append(px(i+1)%mod)
        output.append_line(fixed)
        output.append_line("\nOf course, this has extraneous data, so he removes the extra 2k values to get:")
        snipped = []
        for i in range(n):
            snipped.append(px(i+1)%mod)
        output.append_line(snipped)
    
    

    #Optional prints for if user wanted to see Lagrange, but did not want to see the original data at the beginning
    if show_lagrange and not show_orig:
        output.append_line("\nFor completeness, here is the lagrange interpolation she used to get P(x) originally:\n")
        p = lagrange(point_list,mod,show_lagrange)
        output.append_line("This results in: " + str(p))
    output.append_line("")
    return output.print_output()
