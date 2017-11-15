results = []

def permAux(n, counts, stack, loc):
	for i in range(0,3):
		if counts[i] > 0:
			print(counts)
			stack.append(i)
			#print("Stack after Appending: %s" % stack)
			modified = copier(counts)
			modified[i] -= 1
			if not zeroCheck(modified):
				#print("Recursing, i is %s, stack is %s" % (i, stack))
				tempStack = copier(stack)
				stack = tempStack
				permAux(n, modified, tempStack, loc)
				stack = [loc]
			else:
				print("Appending this stack to results: %s" % stack)
				results.append(stack)

def zeroCheck(arr):
	for i in range(0, len(arr)):
		if arr[i] != 0:
			return False
	return True

def copier(arr):
	result = []
	for i in range(0, len(arr)):
		result.append(arr[i])
	return result

def perm(counts):
	for i in range(0,3):
		if counts[i] > 0:
			stack = []
			stack.append(i)
			modified = copier(counts)
			modified[i] -= 1
			permAux(len(counts), modified, [i], i)


perm([2,1,1])
#print(results)