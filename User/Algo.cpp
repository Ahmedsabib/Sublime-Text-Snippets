// ==================== Binary Exponention ===========================

// iterative method
// (a^b)
long long binpow(long long a, long long b) {
	long long res = 1;
	while (b > 0) {
		if (b&1) {
			res = res * a;
		}
		a = a * a;
		b >>= 1;
	}
	return res;
}


// iterative method
// (a^b)%m
long long binpow(long long a, long long b, long long m) {
	a %= m;
	long long res = 1;
	while (b > 0) {
		if (b&1) {
			res = res * a % m;
		}
		a = a * a % m;
		b >>= 1;
	}
	return res;
}

// ========= GCD (Greatest Common Divisors) && LCM (Least Common Multiple) ================

// Time Complexity O(log min(a, b))
// recursive method
// Method 1
long long gcd(long long a, long long b) {
	return b ? gcd(b, a%b) : a;
}

// LCM
long long lcm(long long a, long long b) {
	return a / gcd(a, b) * b;
}


// ========== Integer Factorization ====================

// Prime Factorization
// Time Complextity O(sqrt(n))
// Method 2 (Most Effecient)
vector<long long> prime_factorization(long long n) {
	vector<long long> factorization;
	while (n%2 == 0) {
		factorization.push_back(2);
		n /= 2;
	}
	for (long long d = 3; d * d <= n; d += 2) {
		while (n%d == 0) {
			factorization.push_back(d);
			n /= d;
		}
	}
	if (n > 1) {
		factorization.push_back(n);
	}
	return factorization;
}

// Finding All Factors of an Integer
// Time Complextity O(sqrt(n))
vector<long long> divisors(long long n) {
	vector<long long> p;
	for (int i = 1; i <= sqrt(n); ++i) {
		if (n%i == 0) {
			if (n/i == i) {
				p.push_back(i);
				continue;
			}
			else {
				p.push_back(n/i);
				p.push_back(i);
			}
		}
	}
	return p;
}

// ==================== Binary Search ============================

// Time Complexity O(log n)
// Method 1
int binary_search(vector<int> array, value) {
	int low = 0, high = (int)array.size() - 1;
	int mid;
	while (low <= high) {
		mid = low + (high - low)/2;
		if (array[mid] == value) {
			return mid;
		}
		if (array[mid] < value) {
			low = mid + 1;
		}
		else {
			high = mid - 1;
		}
	}
	return -1;
}

// lower_bound returns a pointer to the first array element whose value is at least x. (<= x)
auto a = lower_bound(array, array+n, x);

// upper_bound returns a pointer to the first array element whose value is larger than x. (> x)
auto b = upper_bound(array, array+n, x);

// equal_range returns both above pointers.
auto r = equal_range(array, array+n, x);
cout << r.second-r.first << "\n";

// defining int in lower/upper boun(not a pointer)
int k = lower_bound(array, array+n, x) - array;
int l = upper_bound(array, array+n, x) - array;


// =================== Depth First Search (DFS) =======================

// Recursive Method

// Method 1 (for each loop)
// Time Complexity O(n+m)
// Space Complexity O(n)
// n is the number of vertex and
// m is the number of edges
vector<long long> adj[maxn];
bool visited[maxn];
void dfs(long long node) {
	visited[node] = true;
	// using for each loop
	for (auto& edges: adj[node]) {
		if (!visited[edges]) {
			dfs(edges);
		}
	}
}


// Connected Components

// Time Complexity O(n+m)
// n is the number of vertex and
// m is the number of edges
int n;
vector<long long> adj[maxn];
bool visited[maxn];
vector<long long> a;
void dfs(long long node) {
	visited[node] = true;
	a.push_back(node);
	for (int edges = 0; edges < (int)adj[node].size(); ++edges) {
		if (!visited[adj[node][edges]]) {
			dfs(adj[node][edges]);
		}
	}
}
void connected_components() {
	for (int i = 0; i < n; ++i) {
		visited[i] = false;
	}
	for (int i = 0; i < n; ++i) {
		if (!visited[i]) {
			a.clear();
			dfs(i);
			cout << "Component: ";
			for (int j = 0; j < (int)a.size(); ++j) {
				cout << ' ' << a[j];
			}
			cout << '\n';
		}
	}
}

// Connectivity Checking

vector<long long> adj[maxn];
bool visited[maxn];
void dfs(long long node) {
	visited[node] = true;
	// using for loop
	for (int edges = 0; edges < (int)adj[node].size(); ++edges) {
		if (visited[adj[node][edges]] == false) {
			dfs(adj[node][edges]);
		}
	}
}
bool is_connected(long long n) {
	dfs(1);
	for (int i = 1; i <= n; ++i) {
		if (!visited[i]) {
			return false;
		}
	}
	return true;
}
void add_edges(long long u, long long v) {
	adj[u].push_back(v);
	adj[v].push_back(u);
}
int main() {
	int n = 4;
	add_edges(1, 2);
	add_edges(1, 3);
	add_edges(2, 3);
	add_edges(3, 4);
	if (is_connected(n)) {
		cout << "YES" << '\n';
	}
	else {
		cout << "NO" << '\n';
	}
}


// Bipartiteness checking

vector<long long> adj[maxn];
bool visited[maxn];
vector<long long> color(maxn);
void add_edges(long long u, long long v) {
	adj[u].push_back(v);
	adj[v].push_back(u);
}
bool bipartite(long long node) {
	for (auto& edges: adj[node]) {
		// if vertex edges is not explored before
		if (!visited[edges]) {
			// mark present vertex as visited
			visited[edges] = true;
			color[edges] = !color[node];  //mark color opposite to its parents
			if (!bipartite(edges)) {
				return false;
			}
		}
		// if two adjacent are colored with same color then the graph is not bipartite
		else if (color[edges] == color[node]) {
			return false;
		}
	}
	return true;
}
int main() {
	add_edges(3, 2);
	add_edges(1, 4);
	add_edges(2, 1);
	add_edges(5, 3);
	add_edges(6, 2);
	add_edges(3, 1);
	visited[1] = true;
	color[1] = 0;
	if (bipartite(1)) {
		cout << "Graph is Bipartite" << '\n';
	}
	else {
		cout << "Graph is not Bipartite" << '\n';
	}
}

// =========== Shortest Path Graph Algorithm ==========

// Single-source shortest paths

// Bellman-Ford Algorithm
// Not ideal, Dijkstra is better because of it's O((n + m) log m)
// Better than Dijkstra if there are negative weights, edges and cycles

// Time Complexity O(nm)
// Space Complexity O(n)
const int inf = (int)1e9 + 7;
// takes three elements in the edges
vector<vector<long long>>	edge;
void bellman_ford(long long node) {
	// firstly, every distance is infinity
	vector<long long> dist(n, inf);
	// source nodes distance to itself is 0
	dist[node] = 0;
	// checking for the shortest path distance
	for (int i = 0; i < n - 1; ++i) {
		for (auto& x: edge) {
			long long u, v, w;
			u = x[0];
			v = x[1];
			w = x[2];
			dist[v] = min(dist[v], w + dist[u]);
		}
	}
	// checking for negative cycles
	for (int i = 0; i < n - 1; ++i) {
		for (auto& x: edge) {
			long long u, v, w;
			u = x[0];
			v = x[1];
			w = x[2];
			if (dist[u] != inf && dist[u] + w < dist[v]) {
				cout << "negative cycle found" << '\n';
				return;
			}
		}
	}
	// if there is a negative cycle, then the shortest path cannot be found
	// else print the answer
}

// All-pairs shortest paths

// Floyd-Warshall Algorithm
// This technique can be applied if : n <= 500

// For positive weight
// Time Complexity O(n^3)
// Space Complexity O(n^2)
const int inf = (int)1e9 + 7;
long long adj[n][n];
long long a[n][n];
void floyd_warshall() {
  // First, the algorithm initializes distance using 
  // the adjacency matrix adj of the graph
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        adj[i][j] = 0;
      }
      else if (a[i][j]) {
        adj[i][j] = a[i][j];
      }
      else {
        adj[i][j] = inf;
      }
    }
  }
  // Then the shortest distances can be found as follows
  // Remeber the Sequence of looping -> k - i - j
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        adj[i][j] = min(adj[i][j], adj[i][k] + adj[k][j])
      }
    }
  }
}

// =========== Bitmasking ==================

// And Operation
// sets 1 when both x and y are 1
10110 (22)
11010 (26) &

10010 (18)
// x is divisible by 2^k when (x & (2^k - 1) == 0)


// OR Operation
// sets 1 when x or y is 1
10110 (22)
11010 (26) |

11110 (30)


// XOR Operation
// sets 1 when only x or y, one of them is 1
10110 (22)
11010 (26) ^

01100 (12)


// XOR Adding numbers
// The bitwise OR of two numbers is just the sum of those two numbers 
// if there is no carry involved, otherwise you just add their bitwise AND
// Let’s say, we have a=5(101) and b=2(010), since there is no carry involved, 
// their sum is just a|b
// Now, if we change ‘a’ to 6 which is 110 in binary, 
// their sum would change to a|b + a&b since there is a carry involved.
// Time Complexity - O(log b)
// Space Complexity - O(1)
int add(int a, int b) {
  while (b != 0) {
    int carry = a & b;
    a ^= b;
    b = carry << 1;
  }
  return a;
}

// XOR finding the element that occured one time in an array
int find_num(vector<int> a, int n) {
  int res = 0;
  for (int i = 0; i < n; ++i) { // array length n
    res ^= a[i];
  }
  return res;
}


// XOR Finding missing number from an array which contains numbers from 1...n (length - n-1)
// Time Complexity - O(n)
// Space Complexity - O(1)
int get_missing_num(vector<int> a, int n) {
  int x1 = a[0], x2 = 1;
  for (int i = 1; i < n; ++i) {
    x1 ^= a[i];
  }
  for (int i = 2; i <= n+1; ++i) {
    x2 ^= i;
  }
  return x1 ^ x2;
}

// XOR Swapping number
// Time Complexity - O(1)
// Space Complexity - O(1)
int swap_elem(int a, int b) {
  a = a ^ b;
  b = a ^ b;
  a = a ^ b;
  return {x, y};
}

// XOR of all elements from 1 to n
int compute_xor(int n) {
  if (n%4 == 0) {
    return n;
  }
  else if (n%4 == 1) {
    return 1;
  }
  else if (n%4 == 2) {
    return n+1;
  }
  else {
    return 0;
  }
}

// XOR if x is a power of 2
bool is_power_of_two(int x) {
  return x && (!(x & (x - 1)));
}

// XOR to count set bits in a number
int count_set_bits(int x) {
  int cnt = 0;
  while (x) {
    x &= (x - 1);
  }
  return cnt;
}

// Not Operation
// sets every bit inverted
x -> ~x


// Bit Shifts
x << k // appends k zero bits at the right
x >> k // removes k last bits from the right

// Again,
x << k // x * 2^k
x >> k // x / 2^k

// Applications
1 << k // has a 1 bit in the position k and other bits are 0

/**
 * 1. x & (1 << k) -> when kth bit of a number x is 1 (boolean)
 * 2. x | (1 << k) -> sets the kth bit of x to 1
 * 3. x & ~(1 << k) -> sets the kth bit of x to 0
 * 4. x ^ (1 << k) -> inverts the kth bit of x
 * 5. x & (x - 1) -> sets the last bit of x to 0
 * 6. x & ~x -> sets all 1 bits to 0, except for the last one bit
 * 7. x | (x - 1) -> inverts all the bits after the last one bit
 * 8. x & (x - 1) == 0 -> when x is a power of two
 **/

/**
 * __builtin_popcount(x) -> the number of ones in the number
 * __builtin_popcountll(x) -> the number of ones in the number (long long)
 * __builtin_clz(x) -> the number of zeros at the beginning of the number
 * __builtin_ctz(x) -> the number of zeros at the end of the number
 * __builtin_parity(x) -> the parity (even or odd) of the number of ones (if odd then 1, else 0)
**/


// Set Operations
// intersection -> a & b
// union -> a | b
// complement -> ~a
// difference a\b -> a & (~b)


// Following code prints all element belongs to the set
// Following code prints the bit representation of an int x
for (int i = 0; i < 32; ++i) {
	if (x & (1 << i)) {
		cout << i << ' ';
	}
}


// Following code goes through the subsets of (0, 1, .. n-1)
for (int i = 0; i < (1 << n); ++i) {
	// process subset i
}


// Following code goes through the subsets with exactly k elements
for (int i = 0; i < (1 << n); ++i) {
	if (__builtin_popcount(b) == k) {
		// process subset i
	}
}


// Following code goes through the subsets of a set x
int b = 0;
do {
	// process subset b
} while (b = (b - x) & x);
