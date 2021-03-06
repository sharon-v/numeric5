""" assignment 5
    team members: Hadar Amsalem, Sharon Vazana
    git link: https://github.com/sharon-v/numeric5.git"""


# Linear
# Polynomial
# Lagrange
# Neville


def linear(points, XtoFind):
    """
    :param points: get array of points
    :param XtoFind: x we want to find
    :return: Approach in the Linear method
    """
    p1 = points[0]
    for p2 in range(1, len(points)):
        if p1[0] <= XtoFind <= points[p2][0]:
            print("*** Linear_Interpolation ***")
            print("Approximate value = ", makeLine(p1, points[p2], XtoFind))
            print("---------------------------------")
            return
        p1 = points[p2]


def makeLine(point1, point2, XtoFind):
    """
    :param point1: get point one
    :param point2: get point two
    :param XtoFind: x we want to find
    :return: Approximate value between two points
    """
    y1 = point1[1]
    y2 = point2[1]
    x1 = point1[0]
    x2 = point2[0]
    return (((y1 - y2) * XtoFind) / (x1 - x2)) + ((y2 * x1 - y1 * x2) / (x1 - x2))


def polynomial(points, XtoFind):
    """
    :param points: array of points
    :param XtoFind: x we want to find
    :return: Approach in the Polynomial method
    """
    a, b = makePolynomialMat(points)
    # check
    copyA = copyMat(a)
    copyB = copyMat(b)
    copyA, copyB = createDominantDiagonal(copyA, copyB)
    if (copyA is not None) and (copyB is not None):
        a = copyA
        b = copyB
    #
    coefficients = gaussSeidelIter(a, b)
    print("*** Polynomial_Interpolation ***")
    print("Approximate value = ", getCoefficientsCalcY(coefficients, XtoFind))
    print("---------------------------------")


def getCoefficientsCalcY(coefficients, XtoFind):
    """
    :param coefficients: coefficients of poly
    :param XtoFind: x we want to find
    :return: Approximate value in the Polynomial method
    """
    mySum = 0
    for x in range(len(coefficients)):
        mySum += coefficients[x][0] * (XtoFind ** x)
    return mySum


def gaussSeidelIter(a, b):
    """
    :param a: coefficient matrics
    :param b: solution vector
    :return: prints Gauss Seidel method iterations
    """
    epsilon = 0.00001
    iteration = 0
    flag = True
    prevX = makeMatrics(len(a))  # start as zero vector
    currentX = makeMatrics(len(a))
    matA, vectorB = isolateVariables(a, b)
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        prevX[0][0] = currentX[0][0]
        if iteration >= 100:
            print("The system can't converge :(")
            break
        for i in range(len(a)):
            j = 0
            currentX[i][0] = vectorB[i][0]
            while j < len(a[0]):
                if j is not i:
                    currentX[i][0] += matA[i][j] * currentX[j][0]
                j += 1
        iteration += 1
    return currentX


def makeMatrics(row, col=1):
    """
    :param row: get rows of matrix
    :param col: get columns of matrix
    :return: return zero matrix
    """
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


def makePolynomialMat(points):
    """
    :param points: array of points
    :return: Polynomial Mat
    """
    size = len(points)
    newMat = makeMatrics(size, size)
    newB = makeMatrics(size, size)
    for i in range(size):
        xi = points[i][0]
        for j in range(size):
            newMat[i][j] = xi ** j
        newB[i][0] = points[i][1]
    return newMat, newB


def isolateVariables(a, b):
    """
    :param a: coefficient matrics
    :param b: solution vector
    :return: new matrics & vector after isolation of pivot in each row
    """
    vectorB = makeMatrics(len(a))
    matA = makeMatrics(len(a), len(a[0]))
    for i in range(len(a)):
        vectorB[i][0] = b[i][0] / a[i][i]
        j = 0
        while j < len(a[0]):
            if j is i:
                matA[i][i] = 1
            else:
                matA[i][j] -= a[i][j] / a[i][i]
            j += 1

    return matA, vectorB


def lagrange(points, XtoFind):
    """
    :param points: array of points
    :param XtoFind: x we want to find
    :return: Approach in the Lagrange method
    """
    i = 0
    mySum = 0
    while i < len(points):
        xi = points[i][0]
        yi = points[i][1]
        myL = 1
        j = 0
        while j < len(points):
            if j is not i:
                xj = points[j][0]
                myL *= (XtoFind - xj) / (xi - xj)
            j += 1
        mySum += myL * yi
        i += 1

    print("*** Lagrange_Interpolation ***")
    print("Approximate value = ", mySum)
    print("---------------------------------")


def neville(points, XtoFind):
    """
    :param points: array of points
    :param XtoFind: x we want to find
    :return: Approach in the Neville method
    """
    if len(points) < 4:
        print("Can't use Neville's algorithm with less than 4 points...")
        return
    print("*** Neville's_Algorithm ***")
    print("Approximate value = ", recurssiveNeville(points, 0, len(points) - 1, XtoFind))
    print("---------------------------------")


def recurssiveNeville(points, m, n, XtoFind):
    """
    :param points: array of points
    :param m: start of range
    :param n: end of range
    :param XtoFind: x we want to find
    :return: Approximate value in the Neville method
    """
    # stop condition
    if m is n:
        return points[m][1]
    xm = points[m][0]
    xn = points[n][0]
    return (((XtoFind - xm) * recurssiveNeville(points, m + 1, n, XtoFind)) -
            ((XtoFind - xn) * recurssiveNeville(points, m, n - 1, XtoFind))) / (xn - xm)


# dominant diagonal part

def copyMat(A):
    """
    :param A: get matrix
    :return: the same matrix
    """
    B = makeMatrics(len(A), len(A[0]))
    for i in range(len(A)):
        for j in range(len(A[0])):
            B[i][j] = A[i][j]
    return B


def createDominantDiagonal(A, b=None):
    """
    :param A: get matrix
    :param b: get result vector
    :return: matrix with dominant diagonal
    """
    max = 0
    maxIndex = 0
    sum = 0
    for i in range(len(A)):
        for j in range(len(A)):
            sum += abs(A[i][j])
            if abs(A[i][j]) > max:
                max = abs(A[i][j])
                maxIndex = j
        if (sum - max) <= max:
            A = manualSwapCol(A , maxIndex, i)
        else:
            max = 0
            maxIndex = 0
            for j in range(len(A)):
                sum += abs(A[j][i])
                if abs(A[j][i]) > max:
                    max = abs(A[j][i])
                    maxIndex = j
            if rowSum(A[j]) - max <= max:
                A, b = manualSwapRow(A,b, i, maxIndex)
            else:
                print("ERROR - no dominant diagonal")
                return None, None
    return A, b


def manualSwapRow(a, b, r1, r2):
    """
    manaul rows exchange (without e)
    :param a: get matrix
    :param b: get result vector
    :param r1: get number of row to replace
    :param r2: get number of row to replace
    :return: matrix and result vector after replace rows
    """

    if r2 < len(a) and r1 < len(a):
        temp = a[r1]
        a[r1] = a[r2]
        a[r2] = temp
        if b is not None:
            temp = b[r1]
            b[r1] = b[r2]
            b[r2] = temp
    return a, b


def manualSwapCol(a, c1, c2):
    """
    :param a: get matrix
    :param c1: get number of col to replace
    :param c2: get number of col to replace
    :return: matrix after replace cols
    """
    if c2 < len(a) and c1 < len(a):
        for i in range(len(a)):
            temp = a[i][c1]
            a[i][c1] = a[i][c2]
            a[i][c2] = temp
    return a


def rowSum(line):
    """
    :param line: A list od numbers - line for the matrix
    :return: the sum of all the numbers in abs  in the list
    """
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


# end dominant part


def driver():
    points = [[0, 0],
              [1, 0.8415],
              [2, 0.9093],
              [3, 0.1411],
              [4, -0.7568],
              [5, -0.9589],
              [6, -0.2794]]

    XtoFind = 2.5

    # points = [[1, 0.8415],
    #           [2, 0.9093],
    #           [3, 0.1411]]
    #
    # XtoFind = 2.5

    # points = [[1, 1],
    #           [2, 4],
    #           [3, 9]]
    #
    # XtoFind = 3

    # points = [[1, 0],
    #           [1.2, 0.112463],
    #           [1.3, 0.167996],
    #           [1.4, 0.222709]]
    #
    # XtoFind = 1.28

    inp = input("Please Enter: \n0 for Linear Interpolation \n1 for Polynomial Interpolation"
                "\n2 for Lagrange Interpolation \n3 for Neville's Algorithm \nelse for all\n->>>  ")

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Points = ", points)
    print("X = ", XtoFind)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    if inp is '0':
        linear(points, XtoFind)
    elif inp is '1':
        polynomial(points, XtoFind)
    elif inp is '2':
        lagrange(points, XtoFind)
    elif inp is '3':
        neville(points, XtoFind)
    else:
        linear(points, XtoFind)
        polynomial(points, XtoFind)
        lagrange(points, XtoFind)
        neville(points, XtoFind)


driver()
