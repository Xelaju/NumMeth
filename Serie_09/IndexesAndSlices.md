Python indexes and slices for a six-element list.
Indexes enumerate the elements, slices enumerate the spaces between the elements.

												Index Notaion		Slice Notation

Index from rear:    -6  -5  -4  -3  -2  -1      x=[0,1,2,3,4,5]		x[1:]==[1,2,3,4,5]
Index from front:    0   1   2   3   4   5      len(x)==6			x[:5]==[0,1,2,3,4]
                   +---+---+---+---+---+---+    x[0]==0				x[:-2]==[0,1,2,3]
                   | a | b | c | d | e | f |    x[5]==5				x[1:2]==[1]
                   +---+---+---+---+---+---+    x[-1]==5			x[1:-1]==[1,2,3,4]
Slice from front:  :   1   2   3   4   5   :    x[-2]==4			
Slice from rear:   :  -5  -4  -3  -2  -1   :						x[start:end:step]
												y=x[:] 'copy'				Examples:
																	x[::-1]==[5,4,3,2,1,0]
																	x[::-2]==[5,3,1]
																	
Good to Know:
x[a,b,c]

len == length

c  	default is +1. sign of c indicates forward or backward, absolute value of c indicates steps. 
	Default is forward with step size 1. Positive means forward, negative means backward.

a 	when c is positive or blank, default is 0. when c is negative, default is -1.

b 	when c is positive or blank, default is len. when c is negative, default is -(len+1).


In forward direction, starts at 0 and ends at len-1
In backward direction, starts at -1 and ends at -len

IMPORTANT: You will only find elements in this range:

	[-len, -len+1, -len+2, ..., 0, 1, 2,3,4 , len -1]

BUT: The range continues

	[...,-len -2 ,-len-1,-len, -len+1, -len+2, ..., 0, 1, 2,3,4 , len -1, len, len +1, len+2 , ....]

Example:
x=[0,1,2,3,4,5]

x[-2,-6,-1] == x[-2:0:-1] == [4,3,2,1]