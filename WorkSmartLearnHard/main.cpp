#include <bits/stdc++.h>

using namespace std;

#define ii pair<int, int>
#define iii pair< pair<int, int>, int >
#define psi pair<string, int>
#define pis pair<int, string>
#define vi vector<int>
#define msi map<string, int>
#define mci map<char, int>
#define f first 
#define s second	
#define ll long long
#define ull unsigned long long
#define pb push_back
#define mp make_pair
#define sqr(x) x*x
#define PI 3.141592654

bool debug = false;

bool prime[1000001];

//long long Rand(long long l, long long h);
//string RandString(int len);

void SieveOfEratosthenes(int n);
bool isPrime(int n);
bool isPalindrome(string s, int L, int R);
unsigned long long modPow(unsigned long long base, unsigned long long exponent, unsigned long long MOD);
unsigned long long findMMI(unsigned long long n);
void computePrefixHash(string str, int n, unsigned long long int prefix[], unsigned long long power[]);
void computeSuffixHash(string str, int n, unsigned long long int suffix[], unsigned long long power[]);
void findLongestPalindromicString();
// Segment Tree
void update(int k, int l, int r);
void query(int k, int l, int r, int i, int j);

const int oo = 1e7;
const ull MOD=1000000007;
const ll maxn=10111;

ull POW[maxn], prefix[maxn], suffix[maxn];

//ll a[1000001];

ll a[1000001], treeMin[4000001], treeMax[4000001], Min, Max;

string text;
int L[30111], N; //LPS Length Array

ll getHashT(int i,int j) {
    return (prefix[j] - prefix[i - 1] * POW[j - i + 1] + MOD * MOD) % MOD;
}

int bsearch(int a[], int left, int right, float val){
	int mid, pos = -1;
	//if(left>right) return pos;
	while(left<=right){
		mid = (left+right)/2;
		if(a[mid] == val) {
			pos = mid;
			break;
		}
		else if(a[mid] > val){
			right = mid - 1;
		}
		else{
			left = mid+1;
		}
	}
	return pos;
}

std::vector<std::string> taskMaker(std::vector<std::string> source, int challengeId) {
    vector<string> res;
    for(int i=0;i<source.size();i++){
    	size_t found = source[i].find("//DB");
    	if(found != string::npos){
    		string num="";
			for(int j=found+5;source[i][j]!='/';j++){
    			num+=source[i][j];
			}
			int n;
			stringstream ss;
			ss << num;
			ss >> n;
            if(n == challengeId){
	            string tmp = "";
	           	for(int k=0;source[i-1][k]==' ';k++){
	           		tmp+=source[i-1][k];
				}
				tmp+=source[i].substr(found+5+num.size()+2);
				source[i-1] = tmp;
        	}
            source.erase(source.begin()+i);
            i--;
		}
    }
    for(int i=0;i<source.size();i++){
        cout << source[i] << endl;
        res.push_back(source[i]);
    }
    return res;
}

int main(){
	//freopen("input.txt", "r", stdin);
	//freopen("basicstring.out", "w", stdout);
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	//SieveOfEratosthenes(1000000);
	//weekday  = (d+=m<3?y--:y-2,23*m/9+d+4+y/4-y/100+y/400)%7  ;
	cout << "Hello Worldasdfasdfsdfsdf\n";
	return 0;
}

// Random
bool isPrime(int n){
	if(n<2) return false;
	if(n==2) return true;
	if(n%2==0) return false;
	else{
		for(int i=3;i<=sqrt(n);i+=2){
			if(n%i==0) return false;
		}
		return true;
	}
}

// A function to check if a string str is palindrome
// in the range L to R
bool isPalindrome(string str, int L, int R)
{
    // Keep comparing characters while they are same
    while (R > L){
        if (str[L++] != str[R--])
            return false;
        else continue;
	}     
    return true;
}

void SieveOfEratosthenes(int n)
{
    memset(prime, true, sizeof(prime));

    for (int p=2; p*p<=n; p++)
    {
        if (prime[p])
        {
            for (int i=p*2; i<=n; i += p)
                prime[i] = false;
        }
    }

    // Print all prime numbers
    /*bool check = true;
    for (int p=2; p<=n; p++)
       if (prime[p]){
    		v1.pb(p);
	   }*/
}

void update(int k, int l, int r){
	if(l==r){
		treeMin[k] = a[l];
		treeMax[k] = a[l];
		return;
	}
	int m = (l+r)/2;
	update(2*k, l, m);
	update(2*k+1, m+1, r);
	treeMin[k] = min(treeMin[k*2+1], treeMin[k*2]);
	treeMax[k] = max(treeMax[k*2+1], treeMax[k*2]);
}

void query(int k, int l, int r, int i, int j){
	if(l>j || r < i) return;
	if(i<=l && j >= r){
		Min = min(Min, treeMin[k]);
		Max = max(Max, treeMax[k]);
		return;
	}
	int m = (l+r)/2;
	query(2*k, l, m, i, j);
	query(2*k+1, m+1, r, i, j);
}


// A Function to find pow (base, exponent) % MOD
// in log (exponent) time
unsigned long long modPow(unsigned long long base,
                              unsigned long long exponent, unsigned long long MOD)
{
	// Khong de quy
	/*base %= MOD;
    ull result = 1;
    while (exp > 0) {
		if (exp & 1) result = (result * base) % MOD;
    	base = (base * base) % MOD;
    	exp >>= 1;
	}
    return result;*/
    if (exponent == 0)
        return 1;
    if (exponent == 1)
        return base;
 
    unsigned long long int temp = modPow(base, exponent/2, MOD);
 
    if (exponent %2 ==0)
        return (temp%MOD * temp%MOD) % MOD;
    else
        return ((( temp%MOD * temp%MOD)%MOD) * base%MOD) % MOD;
}
 
// A Function to calculate Modulo Multiplicative Inverse of 'n'
unsigned long long findMMI(unsigned long long n)
{
    return modPow(n, MOD-2, MOD);
}

// A Function to calculate the prefix hash
void computePrefixHash(string str, int n, unsigned long long
                       prefix[], unsigned long long power[])
{
    prefix[0] = 0;
    prefix[1] = str[0];
 
    for (int i=2; i<=n; i++)
        prefix[i] = (prefix[i-1]%MOD +
                    (str[i-1]%MOD * power[i-1]%MOD)%MOD)%MOD;
 
    return;
}
 
 
// A Function to calculate the suffix hash
// Suffix hash is nothing but the prefix hash of
// the reversed string
void computeSuffixHash(string str, int n,
                       unsigned long long suffix[],
                       unsigned long long power[])
{
    suffix[0] = 0;
    suffix[1] = str[n-1];
 
    for (int i=n-2, j=2; i>=0 && j<=n; i--,j++)
        suffix[j] = (suffix[j-1]%MOD +
                     (str[i]%MOD * power[j-1]%MOD)%MOD)%MOD;
    return;
}

void findLongestPalindromicString()
{
    N = text.size();
    if(N == 0)
        return;
    N = 2*N + 1; //Position count
    L[0] = 0;
    L[1] = 1;
    int C = 1; //centerPosition 
    int R = 2; //centerRightPosition
    int i = 0; //currentRightPosition
    int iMirror; //currentLeftPosition
    int maxLPSLength = 0;
    int maxLPSCenterPosition = 0;
    int start = -1;
    int end = -1;
    int diff = -1;
     
    //Uncomment it to print LPS Length array
    //printf("%d %d ", L[0], L[1]);
    for (i = 2; i < N; i++) 
    {
        //get currentLeftPosition iMirror for currentRightPosition i
        iMirror  = 2*C-i;
        L[i] = 0;
        diff = R - i;
        //If currentRightPosition i is within centerRightPosition R
        if(diff > 0)
            L[i] = min(L[iMirror], diff);
 
        //Attempt to expand palindrome centered at currentRightPosition i
        //Here for odd positions, we compare characters and 
        //if match then increment LPS Length by ONE
        //If even position, we just increment LPS by ONE without 
        //any character comparison
        int cnt = 0;
        while ( ((i + L[i]) < N && (i - L[i]) > 0) && 
            ( ((i + L[i] + 1) % 2 == 0) || 
            (text[(i + L[i] + 1)/2] == text[(i - L[i] - 1)/2] )))
        {
            L[i]++;
            cnt++;
        }
 		cout << cnt << endl;
        if(L[i] > maxLPSLength)  // Track maxLPSLength
        {
            maxLPSLength = L[i];
            maxLPSCenterPosition = i;
        }
 
        //If palindrome centered at currentRightPosition i 
        //expand beyond centerRightPosition R,
        //adjust centerPosition C based on expanded palindrome.
        if (i + L[i] > R) 
        {
            C = i;
            R = i + L[i];
        }
        //Uncomment it to print LPS Length array
        //printf("%d ", L[i]);
    }
    //printf("\n");
    start = (maxLPSCenterPosition - maxLPSLength)/2;
    end = start + maxLPSLength - 1;
    for(i=start; i<=end; i++)
        printf("%c", text[i]);
    printf("\n");
}

/*long long Rand(long long l, long long h)
{
    return l + ((long long)rand() * (RAND_MAX + 1) * (RAND_MAX + 1) * (RAND_MAX + 1) +
                (long long)rand() * (RAND_MAX + 1) * (RAND_MAX + 1) +
                (long long)rand() * (RAND_MAX + 1) +
                rand()) % (h - l + 1);
}

string RandString(int len){
	string s(len, 42);
	static const char alphanum[] =
        "0123456789";
        //"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        //"abcdefghijklmnopqrstuvwxyz";
    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }
	s[len] = 0;
	return s;
}*/
/*int k = 1;
	srand(time(NULL));
	while(k<=50){
		stringstream ss;
		ss << k;
		string g;
		ss >> g;
		string str = "TEST";
		str += g;
		string m = str;
		str += ".in";
		ifstream ifile(str.c_str());
		//ull n;
		//n = Rand(3, 1000000000000000000);
		//cout << n << endl;
		//ifile >> n;
		string s;
		ifile >> s;
		cout << s << endl;
		m += ".out";
		ofstream ofile(m.c_str());
		if(s=="0"){
			ofile << 1;
			k++;
			continue;
		}
		int n = s[s.size()-1]-'0';
		int choice;
		if(s.size()<2){
			choice = n%4;
		}
		else{
			int tmp = s[s.size()-2] - '0';
			//cout << n << endl << tmp << endl;
			if(n==1 || n==5 || n==9){
				if(tmp&1) choice = 3;
				else choice = 1;
			}
			else if(n==0 || n==4 || n==8){
				if(tmp&1) choice = 2;
				else choice = 0;
			}
			else if(n==3 || n==7){
				if(tmp&1) choice = 1;
				else choice = 3;
			}
			else if(n==2 || n==6){
				if(tmp&1) choice = 0;
				else choice = 2;
			}
		}
		switch(choice){
			case 1: ofile << 8; break;
			case 2: ofile << 4; break;
			case 3: ofile << 2; break;
			case 0: ofile << 6; break;
		}
		k++;
	}*/

