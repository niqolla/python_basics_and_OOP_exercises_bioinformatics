import math

# Exercice_1:
# 1) Write a function that returns a float corresponding to the volume of a sphere:

def get_sphere_volume(radius):
    volume=(4/3)*math.pi*(radius**3)
    return float(volume)

# Exercice_2:
# 2) Write a function that calculates and returns an integer corresponding to the
# factorial of an integer (n):
# a) Using recursivity: recursive_factorial(n)
# b) Without using recursivity: factorial(n)

### a)

def recursive_factorial(n):
    if n == 1:
        return int(n)
    else:
        return int(n*recursive_factorial(n-1))

### b)

def factorial(n):
    i = 1
    result = 1
    while i <= n:
        result*=i
        i+=1
    return int(result)

# Exercice_3:
# 3) Write a function for counting up numbers from 0 to n, showing the count up in the
# screen. If parameter odd is set to True, prints only odd numbers
# a) Using recursivity: recursive_count_up(n, odd)
# b) Without using recursivity: count_up(n,odd)

### a)

def recursive_count_up(n, odd=""):
    # odd is set to false if "True" isn't specified esplicitly 
    if odd==True:
        if n >= 0:
            recursive_count_up(n-1, odd)
            if n%2!=0:
                print(n)
    else:
        if n >= 0:
            recursive_count_up(n-1, odd)
            print(n)

### b)

def count_up(n,odd=""):
    i = 1
    if odd==True:
        while i <= n:
            if i%2!=0:
                print(i)
            i+=1
    else:
        print(0)
        while i <= n:
            print(i)
            i+=1

# Exercice_4:
# 4) Find and solve the bugs in the following function:
# def get_final_price(discount_percentage=10, price):
# """Return the final price after applying the discount percentage """
# return ((price + price) * percentage) / 100

def get_final_price(price, discount_percentage=10):
    """Return the final price after applying the discount percentage """
    return float(price - (price*discount_percentage/100))


# # Examples of use:
# # 1.
# get_sphere_volume(4) # 268.082573106329
# type(get_sphere_volume(4)) # float
# # 2.a
# recursive_factorial(7) # 5040
# type(recursive_factorial(7)) # int
# # 2.b
# factorial(7) # 5040
# type(factorial(7)) # int
# # 3.a
# recursive_count_up(7) # prints 0...7 inclusively; odd=false is set implicitly, can be also specified
# recursive_count_up(7,True) # prints the odd numbers in that range: 1, 3, 5, 7
# recursive_count_up(6) # prints 0...6;
# recursive_count_up(6, True) # prints 1, 3, 5
# # 3.b # comments are the same as in 3.a
# count_up(7)
# count_up(7,True)
# count_up(6)
# count_up(6,True)
# # 4.
# get_final_price(50) # 45.0
# get_final_price(50,50) # 25.0
# get_final_price(432,4) # 414.72
