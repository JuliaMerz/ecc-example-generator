import random
import operator
import collections
import itertools
from itertools import izip
from math import fabs


"""
Practice problem generator with solutions for Reed Solomon codes as presented in CS70, Spring 2013 at UC Berkeley.
Written by Quinn Johnson, tex-based additions by Sunil Srinivasan
"""

#EGCD as presented in the notes. Has extra parameter of "show" to show work, and passes up a list containing tuples of (x, y, d, a, b) for tex to use
def egcd(x, y,show, texlist=[]):
        if y == 0:
                if show and show != 'tex':
                        print "x = " + str(x) + ", y = " + str(y)
                        print "d = " + str(x) + ", a = 1, b = 0"
                texlist.append((x, y, x, 1, 0))
                return (x, 1, 0)
        else:
                if show and show != 'tex':
                        print "x = " + str(x) + ", y = " + str(y)
                (d, a, b) = egcd(y, x % y,show, texlist)
                texlist.append((x, y, d, b, a - (x / y) * b))
                if show and show != 'tex':
                        print "d = " + str(d) + ", a = " + str(b) + ", b = " + str(a - (x / y) * b)
                return ((d , b , a - (x / y) * b))


#Function that calculates the multiplicative inverse of a number by calling egcd. 
#Has parameter of show to show work.
def multinv(num, mod,show):
        endstring = ""
        if show == 'tex':
                endstring = " \\\\"
        if num < 0:
                num = num%mod
        texlist = []
        x = egcd(mod, num,show, texlist)[2]
        texlist.reverse()
        if show == 'tex' and texlist:
                temp = "\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n" + "x & y & d & a & b \\\\ \\hline \n"
                for item in texlist:
                        temp = temp + str(item[0]) + " & " + str(item[1]) + " & " + str(item[2]) + " & " + str(item[3]) + " & " + str(item[4]) + "\\\\ \\hline\n"
                temp = temp + "\\end{tabular}" + endstring
                print temp
        if x < 0:
                if show:
                        print str(x) + " = " + str(mod + x) + " mod " + str(mod) + endstring
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
                        print ""
                        if show == 'tex':
                                print "$\\begin{bmatrix}"
                        for row in matrix:
                                if show == 'tex':
                                        temp = ""
                                        for col in range(0, cols):
                                                temp = temp + str(row[col])
                                                if col == cols-1:
                                                        temp = temp + " \\\\"
                                                else:
                                                        temp = temp + " & "
                                        print temp
                                else:
                                        print row
                        if show == 'tex':
                                print "\\end{bmatrix}$ \\\\\\\\"
                        if show != 'tex':
                                print ""
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

#New function added to output polynomial as tex
    def texstr(poly):
        "Return string formatted as ax^{3} + bx^{2} + cx + d"
        res = []
        for po,co in enumerate(poly.coeffs):
            if co:
                if po==0:
                    po = ''
                elif po==1:
                    po = 'x'
                else:
                    po = 'x^{'+str(po)+'}'
                res.append(str(co)+po)
        if res:
            res.reverse()
            res.append
            return ' + '.join(res)
        else:
            return "0"


#Lagrange Interpolation
#You pass in a list of points in the form [[x1,y1],[x2,y2],...], a mod, and a show parameter
#Returns a Polynomial that matches all points, of degree 1 less than number of points.
#Show parameter is there to allow showing of work

def lagrange(points,mod,show):
        endstring = ""
        if show == 'tex':
                endstring = " \\\\"
        #Construct list of Delta polynomials
        delta_list = []
        #Optional statement for showing work
        if show:
                print "Sample points are: " + str(points).replace('[', '(').replace(']', ')') + endstring
        tempstr = str
        if (show == 'tex'):
                tempstr = Polynomial.texstr
        #Go through each point, creating a delta based on the x value
        #Strings are constructed during the process to allow shown work
        for i in range(len(points)):
                numer = Polynomial(1)
                denom = 1
                numer_string = ""
                denom_string = ""
                for j in range(len(points)):
                        if i != j:
                                numer_string = numer_string + "(" + tempstr(Polynomial((-1*points[j][0]),1)) + ") "
                                denom_string = denom_string + "(" + str(points[i][0]) + " - " + str(points[j][0]) + ") "
                                numer = numer * Polynomial((-1 * points[j][0]),1)
                                denom = denom * (points[i][0] - points[j][0])
                denom = multinv(denom, mod,False)
                delta_list.append(numer * denom)
                #Optional statement for showing work
                if show == 'tex':
                        print "$\\Delta_{" + str(points[i][0]) + "} = " + numer_string + " \\times (" + denom_string + ")^{-1} = " + tempstr(numer * denom) + "$ \\\\"
                elif show:
                        print "Delta " + str(points[i][0]) + ": " + numer_string + " * inverse of " + denom_string + " = " + str(numer * denom)
        val = Polynomial(0)
        if show:
                print "\nWe now add these all together, weighting them with the appropriate Y value:" + endstring
        final_string = ""

        #Combining weighted deltas to get the final polynomial
        if show == 'tex':
                final_string = final_string + "$"
        for i in range(len(points)):
                final_string = final_string + "(" + str(points[i][1]) + ")(" + tempstr(delta_list[i]) + ")"
                if i < len(points) - 1:
                        final_string = final_string + " + "
                val = val + (delta_list[i] * points[i][1])
        if show == 'tex':
                final_string = final_string + "$ \\\\"
        if show:
                print final_string
        return val.mod(mod)


#Solver that takes in a matrix, modulus, and show parameter and solves it.
#For simplicity, makes the assumption that any free variables equal 1
#Returns a vector of coefficients
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
        endstring = ""
        if show == 'tex':
                endstring = " \\\\"
        p_list = []
        #print the form of Q(x) and E(x) if user wants
        if show:
                print "Because Q(x) is of degree " + str(q_deg) + ", and E(x) is of degree " + str(e_deg) + ", he knows:" + endstring
                if show == 'tex':
                        q_temp = q_deg
                        q_string = ""
                        while not q_temp == 0:
                                q_string += "(a_{" + str(q_temp) + "})x^{" + str(q_temp) + "} + "
                                q_temp = q_temp - 1
                        q_string += "(a_0) \\mod " + str(mod)
                        
                        e_temp = e_deg
                        e_string = "x^{" + str(e_temp) + "} + "
                        e_temp = e_temp - 1
                        while not e_temp == 0:
                                e_string += "(b_{" + str(e_temp) + "})x^{" + str(e_temp) + "} + "
                                e_temp = e_temp - 1
                        e_string += "(b_0) \\mod " + str(mod)
                        print "$R(x)$ is defined by the points Bob received." + endstring
                        print "$Q(x) = " + q_string + "$" + endstring
                        print "$E(x) = " + e_string + "$" + endstring + endstring
                else:
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
                        print "R(x) is defined by the points Bob received." + endstring
                        print "Q(x) = " + q_string + endstring
                        print "E(x) = " + e_string + "\n" + endstring

        #For each value of x, calculate coefficients
        for i in range(1,length+1):
                qvec = []
                evec = []
                for j in range(q_deg+1):
                        qvec.append(pow(i,j)%mod)
                for k in range(e_deg+1):
                        evec.append(pow(i,k)%mod)
                
                if show == 'tex':
                        print "Equation " + str(i) + " for $x = " + str(i) + "$:" + endstring
                        num = "(" + str(i) + ")"
                        print "$Q" + num + " = " + q_string.replace('x', num) + "$" + endstring
                        q_poly = Polynomial.texstr(Polynomial(qvec)).replace('x^', 'a_') + "a_0"
                        q_poly = q_poly.replace('x', 'a_1')
                        print "$\\hspace*{26pt} = " + q_poly + " \\mod " + str(mod) + "$" + endstring
                        print "$E" + num + " = " + e_string.replace('x', num) + "$" + endstring
                        e_poly = Polynomial.texstr(Polynomial(evec)).replace('x^{' + str(e_deg) + '}', "")
                        e_poly = e_poly.replace('x^', 'b_') + "b_0"
                        e_poly = e_poly.replace('x', 'b_1')
                        print "$\\hspace*{26pt} = " + e_poly + " \\mod " + str(mod) + "$" + endstring
                        print "$Q" + num + " = R" + num + "E" + num + "$: $" + q_poly + " = " + str(rec[i-1]) + "(" + e_poly + ") \\mod " + str(mod) + "$" + endstring
                elif show:
                        print "Equation " + str(i) + " for x = " + str(i) + ":"
                        num = "(" + str(i) + ")"
                        print "Q" + num + "\t= " + q_string.replace('X', num)
                        q_poly = str(Polynomial(qvec)).replace('X^', 'a') + "a0"
                        q_poly = q_poly.replace('X', 'a1')
                        print "\t= " + q_poly + " mod " + str(mod)
                        print "E" + num + "\t= " + e_string.replace('X', num)
                        e_poly = str(Polynomial(evec)).replace('X^' + str(e_deg), "")
                        e_poly = e_poly.replace('X^', 'b') + "b0"
                        e_poly = e_poly.replace('X', 'b1')
                        print "\t= " + e_poly + " mod " + str(mod)
                        print "Q" + num + "=R" + num + "E" + num + ":\t" + q_poly + " = " + str(rec[i-1]) + "(" + e_poly + ") mod " + str(mod)

                for l in range(len(evec)):
                        evec[l] = (evec[l] * rec[i-1])%mod
                
                if show == 'tex':
                        e_poly = Polynomial.texstr(Polynomial(evec)).replace('x^{' + str(e_deg) + '}', "")
                        e_poly = e_poly.replace('x^', 'b_') + "b_0"
                        e_poly = e_poly.replace('x', 'b_1')
                        print "\\hspace*{26pt}$\\rightarrow" + q_poly + " = " + e_poly + " \\mod " + str(mod) + "$" + endstring
                elif show:
                        e_poly = str(Polynomial(evec)).replace('X^' + str(e_deg), "")
                        e_poly = e_poly.replace('X^', 'b') + "b0"
                        e_poly = e_poly.replace('X', 'b1')
                        print "\t=>\t" + q_poly + " = " + e_poly + " mod " + str(mod)

                qvec.reverse()
                evec.reverse()
                #Combine into one vector
                vec = qvec
                final_eq = ""
                if show == 'tex':
                        q_temp = q_deg
                        for x in vec:
                                final_eq += str(x) + "a_{" + str(q_temp) + "} + "
                                q_temp = q_temp - 1
                        e_temp = e_deg - 1
                elif show:
                        q_temp = q_deg
                        for x in vec:
                                final_eq += str(x) + "a" + str(q_temp) + " + "
                                q_temp = q_temp - 1
                        e_temp = e_deg - 1

                for j in range(1,e_deg+1):
                        app_val = ((-1)*evec[j])%mod
                        vec.append(app_val)
                        if show == 'tex':
                                final_eq += str(app_val) + "b_{" + str(e_temp) + "}"
                                if not e_temp == 0:
                                    final_eq += " + "
                                e_temp = e_temp - 1
                        elif show:
                                final_eq += str(app_val) + "b" + str(e_temp)
                                if not e_temp == 0:
                                    final_eq += " + "
                                e_temp = e_temp - 1

                vec.append(evec[0])
                
                if show == 'tex':
                        final_eq += " = " + str(evec[0]) + " \\mod " + str(mod)
                        print "\\hspace*{26pt}$\\rightarrow" + final_eq + "$" + endstring
                        pause("\\\\")
                elif show:
                        final_eq += " = " + str(evec[0]) + " mod " + str(mod)
                        print "\t=>\t" + final_eq + "\n\n"
                        pause()

                p_list.append(vec)
        return p_list

#Gets intial values from the user
def get_user_input():
        try:
                print "Please enter variables relevant to the Reed-Solomon problem:"

                show_inv = ""
                show_rref = ""
                show_lagrange = ""
                show_orig = ""
                show_bw = ""

                n = int(raw_input("Length of original message (n in notes): \n"))
                k = int(raw_input("Number of errors to protect against (k in notes): \n"))
                mod = int(raw_input("Prime number we are working modulo (must be prime and bigger than n + 2k)\nEnter \"0\" to default to next highest prime: \n"))
                if mod == 0:
                        i = n + 2*k + 1
                        while not prime(i):
                                i += 1
                        mod = i
                        print "Calculations will be done modulo " + str(mod) + "\n"
                
                print "\nNow enter what steps you want to see covered in detail. All set to n provides a minimal solution."

                #Demand y, Y, tex, n, or N
                while (show_inv != 'y' and show_inv != 'Y' and show_inv != 'tex' and show_inv != 'n' and show_inv != 'N'):
                        show_inv = str(raw_input("Multiplicative inverse / EGCD? (y/tex/n): \n"))
                if show_inv == 'n' or show_inv == 'N':
                        show_inv = False

                while (show_rref != 'y' and show_rref != 'Y' and show_rref != 'tex' and show_rref != 'n' and show_rref != 'N'):
                        show_rref = str(raw_input("Row Reduction to solve system of equations? (y/tex/n): \n"))
                if show_rref == 'n' or show_rref == 'N':
                        show_rref = False

                while (show_lagrange != 'y' and show_lagrange != 'Y' and show_lagrange != 'tex' and show_lagrange != 'n' and show_lagrange != 'N'):
                        show_lagrange = str(raw_input("Lagrange Polynomial Interpolation? (y/tex/n): \n"))
                if show_lagrange == 'n' or show_lagrange == 'N':
                        show_lagrange = False

                while (show_orig != 'y' and show_orig != 'Y' and show_orig != 'tex' and show_orig != 'n' and show_orig != 'N'):
                        show_orig = str(raw_input("Do you want to see the original polynomial/message at the beginning? (y/tex/n): \n"))
                if show_orig == 'n' or show_orig == 'N':
                        show_orig = False

                while (show_bw != 'y' and show_bw != 'Y' and show_bw != 'tex' and show_bw != 'n' and show_bw != 'N'):
                        show_bw = str(raw_input("Demonstrate the Berlekamp-Welch algorithm to construct a system of equations? (y/tex/n): \n"))
                if show_bw == 'n' or show_bw == 'N':
                        show_bw = False
                
                print "\n\n"
                
                
        except:
                print "Input error"
                exit()
        return [n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw]


#Pauses for user input to continue
def pause(endstring = "Press Enter to continue.."):
        raw_input(endstring)


#Main interactive function
def run_interactive():
        random.seed()
        
        #Get values from user
        [n, k, mod, show_inv, show_rref, show_lagrange, show_orig, show_bw] = get_user_input()
        length = n + (k*2)
        orig = [0]*n
        rec = [0]*length
        point_list = []
        #tex-compatibility lines
        endstring = ""
        pausestring = "Press Enter to continue...\n"
        if show_inv == 'tex' and show_rref == 'tex' and show_lagrange == 'tex' and ((not show_orig) or (show_orig == 'tex')):
                endstring = "\\\\"
                pausestring = "\\\\"
        #Generate original message
        for i in range(n):
                orig[i] = random.randint(0,mod-1)
                point_list.append([i+1,orig[i]])

        print "For this problem, Alice wants to send a message of length " + str(n) + " to Bob." + endstring
        print "She also wants to ensure that it can survive up to " + str(k) + " general errors." + endstring
        print "She realizes that through the Reed-Solomon scheme, she can do so with a message of length " + str(length) + " over GF(" + str(mod) +")." + endstring
        if (show_orig):
                print "\nHer message is:\n " + str(orig).replace('[', '(').replace(']', ')') + "\n" + endstring
        
        #Calculating degrees of polynomials
        q_deg = n + k -1
        e_deg = k
        p_deg = n - 1

        print "We know that if we create an original polynomial P(x) of degree " + str(p_deg) + " passing through the points of her message, we can extend it." + endstring
        print "From the extended (and possibly corrupted) message, we can obtain P(x) from two polynomials:" + endstring
        print "The error polynomial, E(x), of degree " + str(e_deg) + " and P(x)E(x) = Q(x), of degree " + str(q_deg) + "." + endstring

        #Many optional print statements if user wanted to see original information before being sent
        if (show_orig):
                print "\nHere is the interpolation she does to find her original polynomial, P(x)\n" + endstring
        
        show_lagrange_early = show_lagrange and show_orig
        if (show_lagrange_early and show_lagrange == 'tex'):
                show_lagrange_early = 'tex'

        #Calculate the original polynomial
        orig_poly = lagrange(point_list,mod,show_lagrange_early)
        if (show_orig == 'tex'):
                print "\nHer original polynomial, obtained through interpolation, is: \n$" + orig_poly.texstr() + "$" + endstring
        elif show_orig:
                print "\nHer original polynomial, obtained through interpolation, is: \n" + str(orig_poly)

        #Creates the message Alice sends
        for i in range(length):
                rec[i] = orig_poly(i+1)%mod

        if (show_orig):
                print "Using this, she sends out the full message: " + str(rec).replace('[', '(').replace(']', ')') + "\n" + endstring
        #Randomizes the message to create errors
        for _ in range(k):
                a = random.randint(0,length-1)
                rec[a] = random.randint(0,mod-1)

        pause(pausestring)

        #Set of print statements and computations if user wants to see Multiplicative Inverse calculations
        if show_inv:
                print "Because you, the user, asked to see the computation of the inverses, we will do so for all values less than our modulus.\n" + endstring
                inv_list = [-1]*mod
                for i in range(1,mod):
                        if inv_list[i] == -1:
                                print "We have not yet calculated the inverse of " + str(i) + " so we will do that with egcd:\n" + endstring
                                a = multinv(i, mod, show_inv)
                                inv_list[a] = i
                                inv_list[i] = a
                                print "From this, we know that " + str(a) + " is " + str(i) + "\'s inverse, and vice-versa.\n" + endstring
                                pause(pausestring)
                        else:
                                print "We have already calculated the inverse of " + str(i) + " so we will skip it. \n" + endstring
                                pause(pausestring)

        #Explain how to construct system of equations:
        #Could possibly be done with the """ """ construct, but that may interfere with doc-strings if I add them.
        print "Bob receives the message: " + str(rec).replace('[', '(').replace(']', ')') + ", which he knows may be wrong in up to " + str(k) + " places.\n" + endstring
        print "With this, he seeks to create a system of equations of the form Q(x) = R(x)E(x) for each point he received." + endstring
        print "He can reconstruct the original P(x) by solving the system." + endstring

        if show_bw:
                print "He knows that Q(x) has " + str(q_deg + 1) + " degrees of freedom, and E(x) has " + str(e_deg) + " degrees of freedom (because the leading coefficient is always 1).\n" + endstring
                print "To do so, he evaluates Q(x) = R(x)E(x) for each value of x he has received." + endstring
                print "From his received message, R(x) corresponds to a received value. For instance, R(1) = " + str(rec[0]) + "." + endstring
                print "We can construct Q(x) and E(x) as follows:\n" + endstring

                if show_bw == 'tex':
                        print "$Q(x) = a_{(n+k-1)}x^{(n+k-1)} + ... + a_1x + a_0$" + endstring
                        print "$E(x) = x^k + b_{(k-1)}x^{(k-1)} + ... + b_1x + b_0$" + endstring
                else:
                        print "Q(x) = a_(n+k-1) * x^(n+k-1) + ... + a_1 * x + a_0" + endstring #unnecessary endstrings but might as well
                        print "E(x) = x^k + b_(k-1) * x^(k-1) + ... + b_1 * x + b_0" + endstring
                print "Note here that the leading coefficient is always 1." + endstring
                print "Evaluating at each point, he constructs the following equations:" + endstring
                pause(pausestring)
        poly_list = bw_alg(length, rec, q_deg, e_deg, mod, show_bw)

        #Display system of equations
        print "By plugging in values for each x he was given, he creates this system of equations." + endstring
        print "(Those without Math54/equivalent experience please look up row-reduction of matrices for solving sets of linear equations, it's a nice way of doing so)." + endstring
        print "Because it is a system of " + str(q_deg + e_deg + 1) + " linear equations in " + str(q_deg + e_deg + 1) + " variables, he represents it as the following matrix:" + endstring
        if show_rref == 'tex':
                print "$\\begin{bmatrix}"
                for row in poly_list:
                        temp = ""
                        for num in enumerate(row):
                                temp = temp + str(num[1])
                                if num[0] != length:
                                        temp = temp + " & "
                                else:
                                        temp = temp + " \\\\"
                        print temp
                print "\\end{bmatrix}$" + endstring
        else:
                for row in poly_list:
                        print row
        
        pause(pausestring)
        print "\nHe then row-reduces this matrix to solve for each unknown: " + endstring

        #Solution taken is one where all free variables are 1
        sols = matrix_solver(poly_list,mod,show_rref)
        print "\nAfter row reduction, he arrives at the following set of values (note that he takes all free variables to be 1): " + endstring
        print str(sols).replace('[', '(').replace(']', ')') + "\n" + endstring

        #Constructing and displaying solutions as 2 sets of coefficients
        qvec = sols[:q_deg+1]
        print "The first " + str(q_deg+1) + " values correspond to the coefficients of Q(x), from highest degree to lowest: \n" + endstring +  str(qvec).replace('[', '(').replace(']', ')') + "\n" + endstring
        qvec.reverse()
        evec = sols[q_deg+1:]
        print "The remaining values correspond to the coefficients of E(x), from highest degree to lowest: \n" + endstring + str(evec).replace('[', '(').replace(']', ')') + "\n" + endstring
        print "Note that for the coefficients of E(x), it is understood that the leading coefficient is 1, and is not in the above list.\n" + endstring
        evec.reverse()
        evec.append(1)
        
        #Creating polynomials
        qx = Polynomial(qvec)
        ex = Polynomial(evec)
        
        pause(pausestring)
        print "From this, he can reconstruct Q(x) and E(x): " + endstring
        if show_orig == 'tex':
                print "$Q(x) = " + qx.texstr() + "$" + endstring
                print "$E(x) = " + ex.texstr() + "$" + endstring
        else:
                print "Q(x) = " + str(qx)
                print "E(x) = " + str(ex)
        
        #Calculating original polynomial
        pause(pausestring)
        px = (qx.div(ex,mod)).mod(mod)
        if show_orig == 'tex':
                print "\nBob then divides $Q(x)$ by $E(x)$ to finally get back to the original polynomial:" + endstring + "\n$P(x) = " + px.texstr() + "$" + endstring
        else:
                print "\nBob then divides Q(x) by E(x) to finally get back to the original polynomial P(x):\n" + str(px) 
        print "\nBecause the original message lies along the values of P(x), Bob now has sufficient information to find any errors in his received copy." + endstring

        #Finding roots of E(x)
        print "Bob also knows that the roots of E(x) correspond to the locations of errors, so he solves:" + endstring
        if show_orig == 'tex':
                print "$" + ex.texstr() + " = 0$" + endstring
        else:
                print str(ex) + " = 0"
        print "and finds that the message was corrupted at positions:"

        roots = []
        for i in range(0,mod):
                val = ex(i)%mod
                if val == 0:
                        roots.append(i)
        
        #If no errors, special message
        if len(roots) == 0:
                print "...Interesting, there were no errors. Lucky Bob." + endstring
        #Otherwise, print out original
        else:
                pause(pausestring)
                err_locs = ""
                for err in roots:
                        err_locs += str(err) + "; "
                print err_locs + "\n" + endstring
                print "So, he recalculates to find that the actual values are:\n" + endstring
                for err in roots:
                        print "Position " + str(err) + ": " + str(px(err)%mod) + ", not " + str(rec[(err-1)%(length)]) + endstring #something strange was going on here with rec[err-1] occasionally, so mod it just to be... safe?
            
                print "\nFinally, he fixes the errors he found to arrive back at the original transmitted message:" + endstring
                fixed = []
                for i in range(length):
                        fixed.append(px(i+1)%mod)
                print str(rec).replace('[', '(').replace(']', ')') + " becomes " + str(fixed).replace('[', '(').replace(']', ')') + endstring
                print "\nOf course, this has extraneous data, so he removes the extra 2k values to get:" + endstring
                snipped = []
                for i in range(n):
                        snipped.append(px(i+1)%mod)
                print str(snipped).replace('[', '(').replace(']', ')') + endstring


        #Optional prints for if user wanted to see Lagrange, but did not want to see the original data at the beginning
        if show_lagrange and not show_orig:
                print "\nFor completeness, here is the lagrange interpolation she used to get P(x) originally:\n" + endstring
                p = lagrange(point_list,mod,show_lagrange)
                if show_lagrange == 'tex':
                        print "This results in: $" + Polynomial.texstr(p) + "$" + endstring
                else:
                        print "This results in: " + str(p)
        print ""

#Call to main function
run_interactive()
