#include <bits/stdc++.h>
using namespace std;
  
// structure for Binary Node
struct Node {
    int key;
    struct Node *right, *left;
};

Node* newNode(int num)
{
    Node* temp = new Node;
    temp->key = num;
    temp->left = NULL;
    temp->right = NULL;
    return temp;
}
  
// To create a Tree with n levels. We always
// insert new node to left if it is less than
// previous value.
Node* createNLevelTree(int arr[], int n)
{
    Node* root = newNode(arr[0]);
    Node* temp = root;
    for (int i = 1; i < n; i++) {
        if (temp->key > arr[i]) {
            temp->left = newNode(arr[i]);
            temp = temp->left;
        }
        else {
            temp->right = newNode(arr[i]);
            temp = temp->right;
        }
    }
    return root;
}
  
// Please refer below post for details of this
// function.
// https:// www.geeksforgeeks.org/a-program-to-check-if-a-binary-tree-is-bst-or-not/
bool isBST(Node* root, int min, int max)
{
    if (root == NULL)
        return true;
  
    if (root->key < min || root->key > max)
        return false;
  
    // Allow only distinct values
    return (isBST(root->left, min,
                  (root->key) - 1)
            && isBST(root->right,
                     (root->key) + 1, max));
}
  
// Returns tree if given array of size n can
// represent a BST of n levels.
bool canRepresentNLevelBST(int arr[], int n)
{
    Node* root = createNLevelTree(arr, n);
    return isBST(root, INT_MIN, INT_MAX);
}
 



//Check if a given Binary Tree is Heap


/* This function counts the number of nodes in a binary tree */
unsigned int countNodes(struct Node* root)
{
    if (root == NULL)
        return (0);
    return (1 + countNodes(root->left) + countNodes(root->right));
}
  
/* This function checks if the binary tree is complete or not */
bool isCompleteUtil (struct Node* root, unsigned int index,
                     unsigned int number_nodes)
{
    // An empty tree is complete
    if (root == NULL)
        return (true);
  
    // If index assigned to current node is more than
    // number of nodes in tree, then tree is not complete
    if (index >= number_nodes)
        return (false);
  
    // Recur for left and right subtrees
    return (isCompleteUtil(root->left, 2*index + 1, number_nodes) &&
            isCompleteUtil(root->right, 2*index + 2, number_nodes));
}
  
// This Function checks the heap property in the tree.
bool isHeapUtil(struct Node* root)
{
    //  Base case : single node satisfies property
    if (root->left == NULL && root->right == NULL)
        return (true);
  
    //  node will be in second last level
    if (root->right == NULL)
    {
        //  check heap property at Node
        //  No recursive call , because no need to check last level
         return (root->key >= root->left->key);
    }
    else
    {
        //  Check heap property at Node and
        //  Recursive check heap property at left and right subtree
        if (root->key >= root->left->key &&
            root->key >= root->right->key)
            return ((isHeapUtil(root->left)) &&
                    (isHeapUtil(root->right)));
        else
            return (false);
    }
}
  
//  Function to check binary tree is a Heap or Not.
bool isHeap(struct Node* root)
{
    // These two are used in isCompleteUtil()
    unsigned int node_count = countNodes(root);
    unsigned int index = 0;
  
    if (isCompleteUtil(root, index, node_count) && isHeapUtil(root))
        return true;
    return false;
}







// A recursive function to construct BST from pre[]. preIndex is used
// to keep track of index in pre[].
struct node* constructTreeUtil( int pre[], int* preIndex, int key,
                                int min, int max, int size )
{
    // Base case
    if( *preIndex >= size )
        return NULL;
   
    struct node* root = NULL;
   
    // If current element of pre[] is in range, then
    // only it is part of current subtree
    if( key > min && key < max )
    {
        // Allocate memory for root of this subtree and increment *preIndex
        root = newNode ( key );
        *preIndex = *preIndex + 1;
          
        if (*preIndex < size)
        {
            // Contruct the subtree under root
            // All nodes which are in range {min .. key} will go in left
            // subtree, and first such node will be root of left subtree.
            root->left = constructTreeUtil( pre, preIndex, pre[*preIndex],
                                        min, key, size );

            // All nodes which are in range {key..max} will go in right
            // subtree, and first such node will be root of right subtree.
            root->right = constructTreeUtil( pre, preIndex, pre[*preIndex],
                                         key, max, size );
        }
    }
   
    return root;
}
  
// The main function to construct BST from given preorder traversal.
// This function mainly uses constructTreeUtil()
struct node *constructTree (int pre[], int size)
{
    int preIndex = 0;
    return constructTreeUtil ( pre, &preIndex, pre[0], INT_MIN, INT_MAX, size );
}
  












//Construct a Binary Search Tree from given postorder
struct node* constructTreeUtil(int post[], int* postIndex,
                         int key, int min, int max, int size)
{
    // Base case
    if (*postIndex < 0)
        return NULL;
     struct node* root = NULL;
  
    // If current element of post[] is in range, then
    // only it is part of current subtree
    if (key > min && key < max)
    {
        // Allocate memory for root of this subtree and decrement
        // *postIndex
        root = newNode(key);
        *postIndex = *postIndex - 1;
  
        if (*postIndex >= 0)
        {
  
          // All nodes which are in range {key..max} will go in right
          // subtree, and first such node will be root of right subtree.
          root->right = constructTreeUtil(post, postIndex, post[*postIndex],
                                          key, max, size );
  
          // Contruct the subtree under root
          // All nodes which are in range {min .. key} will go in left
          // subtree, and first such node will be root of left subtree.
          root->left = constructTreeUtil(post, postIndex, post[*postIndex],
                                         min, key, size );
        }
    }
    return root;
}
  
// The main function to construct BST from given postorder
// traversal. This function mainly uses constructTreeUtil()
struct node *constructTree (int post[], int size)
{
    int postIndex = size-1;
    return constructTreeUtil(post, &postIndex, post[postIndex],
                             INT_MIN, INT_MAX, size);
}
  















Query of type 1 :
Find the range sum on segment tree for output query where range is exit time and entry
time of the rooted tree node. Deduce that the answer is always twice the expected answer
because each node is added twice in segment tree. So reduce the answer by half.

Query of type 2 :
For update query, update the leaf node of segment tree at the entry time and exit time of
the rooted tree node.
Below is the implementation of above approach :
// C++ program for implementation of
// Euler Tour | Subtree Sum.
#include <bits/stdc++.h>
using namespace std;
  
vector<int> v[1001];
vector<int> s;
int seg[1001] = { 0 };

// Value/Weight of each node of tree,
// value of 0th(no such node) node is 0.
int ar[] = { 0, 1, 2, 3, 4, 5, 6 };
  
int vertices = 6;
int edges = 5;
  
// A recursive function that constructs
// Segment Tree for array ar[] = { }.
// 'pos' is index of current node
// in segment tree seg[].
int segment(int low, int high, int pos)
{
    if (high == low) {
        seg[pos] = ar[s[low]];
    }
    else {
        int mid = (low + high) / 2;
        segment(low, mid, 2 * pos);
        segment(mid + 1, high, 2 * pos + 1);
        seg[pos] = seg[2 * pos] + seg[2 * pos + 1];
    }
}
  
/* Return sum of elements in range
   from index l to r . It uses the 
   seg[] array created using segment()
   function. 'pos' is index of current
   node in segment tree seg[].
*/
int query(int node, int start,
          int end, int l, int r)
{
    if (r < start || end < l) {
        return 0;
    }
  
    if (l <= start && end <= r) {
        return seg[node];
    }
  
    int mid = (start + end) / 2;
    int p1 = query(2 * node, start,
                   mid, l, r);
    int p2 = query(2 * node + 1, mid + 1,
                   end, l, r);
  
    return (p1 + p2);

    }
  
/* A recursive function to update the
   nodes which have the given index in
   their range. The following are
   parameters pos --> index of current 
   node in segment tree seg[]. idx -->
   index of the element to be updated.
   This index is in input array. 
   val --> Value to be change at node idx
*/
int update(int pos, int low, int high,
           int idx, int val)
{
    if (low == high) {
        seg[pos] = val;
    }
    else {
        int mid = (low + high) / 2;
  
        if (low <= idx && idx <= mid) {
            update(2 * pos, low, mid,
                   idx, val);
        }
        else {
            update(2 * pos + 1, mid + 1,
                   high, idx, val);
        }
  
        seg[pos] = seg[2 * pos] + seg[2 * pos + 1];
    }
}
  
/* A recursive function to form array
    ar[] from a directed tree .
*/
int dfs(int root)
{
    // pushing each node in vector s
    s.push_back(root);
    if (v[root].size() == 0)
        return root;
  
    for (int i = 0; i < v[root].size(); i++) {
        int temp = dfs(v[root][i]);
        s.push_back(temp);
    }
    return root;


    }
  
// Driver program to test above functions
int main()
{
    // Edges between the nodes
    v[1].push_back(2);
    v[1].push_back(3);
    v[2].push_back(6);
    v[2].push_back(5);
    v[3].push_back(4);
  
    // Calling dfs function.
    int temp = dfs(1);
    s.push_back(temp);
  
    // Storing entry time and exit
    // time of each node
    vector<pair<int, int> > p;
  
    for (int i = 0; i <= vertices; i++)
        p.push_back(make_pair(0, 0));
  
    for (int i = 0; i < s.size(); i++) {
        if (p[s[i]].first == 0)
            p[s[i]].first = i + 1;
        else
            p[s[i]].second = i + 1;
    }
  
    // Build segment tree from array ar[].
    segment(0, s.size() - 1, 1);
  
    // query of type 1 return the
    // sum of subtree at node 1.
    int node = 1;
    int e = p[node].first;
    int f = p[node].second;
  
    int ans = query(1, 1, s.size(), e, f);
  
    // print the sum of subtree
    cout << "Subtree sum of node " << node << " is : " << (ans / 2) << endl;
  
    // query of type 2 return update
    // the subtree at node 6.
    int val = 10;
    node = 6;
     e = p[node].first;
    f = p[node].second;
    update(1, 1, s.size(), e, val);
    update(1, 1, s.size(), f, val);
  
    // query of type 1 return the
    // sum of subtree at node 2.
    node = 2;
  
    e = p[node].first;
    f = p[node].second;
  
    ans = query(1, 1, s.size(), e, f);
  
    // print the sum of subtree
    cout << "Subtree sum of node " << node << " is : " << (ans / 2) << endl;
  
    return 0;
}



Output:
Subtree sum of node 1 is : 21
Subtree sum of node 2 is : 17
Time Complexity : O(q*log(n))





// This function clones a given linked list
// in O(1) space
Node* clone(Node *start)
{
    Node* curr = start, *temp;
  
    // insert additional node after
     // every node of original list
    while (curr)
    {
        temp = curr->next;
  
        // Inserting node
        curr->next = new Node(curr->data);
        curr->next->next = temp;
        curr = temp;
    }
  
    curr = start;
  
    // adjust the random pointers of the
    // newly added nodes
    while (curr)
    {
        curr->next->random = curr->random->next;
  
        // move to the next newly added node by
        // skipping an original node
        curr = curr->next?curr->next->next:curr->next;
    }
  
    Node* original = start, *copy = start->next;
  
    // save the start of copied linked list
    temp = copy;
  
    // now separate the original list and copied list
    while (original && copy)
    {
        original->next =
         original->next? original->next->next : original->next;
  
        copy->next = copy->next?copy->next->next:copy->next;
        original = original->next;
        copy = copy->next;
    }
  
    return temp;
}
  
// Driver code
int main()
{
    Node* start = new Node(1);
    start->next = new Node(2);







    #include <bits/stdc++.h> 
  
int _mergeSort(int arr[], int temp[], int left, int right); 
int merge(int arr[], int temp[], int left, int mid, int right); 
  
/* This function sorts the input array and returns the 
   number of inversions in the array */
int mergeSort(int arr[], int array_size) 
{ 
    int* temp = (int*)malloc(sizeof(int) * array_size); 
    return _mergeSort(arr, temp, 0, array_size - 1); 
} 
  
/* An auxiliary recursive function that sorts the input array and 
  returns the number of inversions in the array. */
int _mergeSort(int arr[], int temp[], int left, int right) 
{ 
    int mid, inv_count = 0; 
    if (right > left) { 
        /* Divide the array into two parts and call _mergeSortAndCountInv() 
       for each of the parts */
        mid = (right + left) / 2; 
  
        /* Inversion count will be sum of inversions in left-part, right-part 
      and number of inversions in merging */
        inv_count = _mergeSort(arr, temp, left, mid); 
        inv_count += _mergeSort(arr, temp, mid + 1, right); 
  
        /*Merge the two parts*/
        inv_count += merge(arr, temp, left, mid + 1, right); 
    } 
    return inv_count; 
} 
  
/* This funt merges two sorted arrays and returns inversion count in 
   the arrays.*/
int merge(int arr[], int temp[], int left, int mid, int right) 
{ 
    int i, j, k; 
    int inv_count = 0; 
  
    i = left; /* i is index for left subarray*/
    j = mid; /* j is index for right subarray*/
    k = left; /* k is index for resultant merged subarray*/
    while ((i <= mid - 1) && (j <= right)) { 
        if (arr[i] <= arr[j]) { 
            temp[k++] = arr[i++]; 
        } 
        else { 
            temp[k++] = arr[j++]; 
  
            /*this is tricky -- see above explanation/diagram for merge()*/
            inv_count = inv_count + (mid - i); 
        } 
    } 
  
    /* Copy the remaining elements of left subarray 
   (if there are any) to temp*/
    while (i <= mid - 1) 
        temp[k++] = arr[i++]; 
  
    /* Copy the remaining elements of right subarray 
   (if there are any) to temp*/
    while (j <= right) 
        temp[k++] = arr[j++]; 
  
    /*Copy back the merged elements to original array*/
    for (i = left; i <= right; i++) 
        arr[i] = temp[i]; 
  
    return inv_count; 
} 
  
/* Driver program to test above functions */
int main(int argv, char** args) 
{ 
    int arr[] = { 1, 20, 6, 4, 5 }; 
    printf(" Number of inversions are %d \n", mergeSort(arr, 5)); 
    getchar(); 
    return 0; 
} 




// We can use stl container list as a double 
// ended queue to store the cache keys, with 
// the descending time of reference from front 
// to back and a set container to check presence 
// of a key. But to fetch the address of the key 
// in the list using find(), it takes O(N) time. 
// This can be optimized by storing a reference 
//     (iterator) to each key in a hash map. 
#include <bits/stdc++.h> 
using namespace std; 
  
class LRUCache { 
    // store keys of cache 
    list<int> dq; 
  
    // store references of key in cache 
    unordered_map<int, list<int>::iterator> ma; 
    int csize; // maximum capacity of cache 
  
public: 
    LRUCache(int); 
    void refer(int); 
    void display(); 
}; 
  
// Declare the size 
LRUCache::LRUCache(int n) 
{ 
    csize = n; 
} 
  
// Refers key x with in the LRU cache 
void LRUCache::refer(int x) 
{ 
    // not present in cache 
    if (ma.find(x) == ma.end()) { 
        // cache is full 
        if (dq.size() == csize) { 
            // delete least recently used element 
            int last = dq.back(); 
  
            // Pops the last elmeent 
            dq.pop_back(); 
  
            // Erase the last 
            ma.erase(last); 
        } 
    } 
  
    // present in cache 
    else
        dq.erase(ma[x]); 
  
    // update reference 
    dq.push_front(x); 
    ma[x] = dq.begin(); 
} 
  
// Function to display contents of cache 
void LRUCache::display() 
{ 
  
    // Iterate in the deque and print 
    // all the elements in it 
    for (auto it = dq.begin(); it != dq.end(); 
         it++) 
        cout << (*it) << " "; 
  
    cout << endl; 
} 
  
// Driver Code 
int main() 
{ 
    LRUCache ca(4); 
  
    ca.refer(1); 
    ca.refer(2); 
    ca.refer(3); 
    ca.refer(1); 
    ca.refer(4); 
    ca.refer(5); 
    ca.display(); 
  
    return 0; 
} 


// mirror of a tree

void mirror(struct node* node) 
{
    if (node) 
    {
        struct node* temp;
        /* do the subtrees */
        mirror(node->left);
        mirror(node->right);
        /* swap the pointers in this node */
        temp  = node->left;
        node->left  = node->right;
        node->right = temp;
     }
}


// Maximum Width of Binary Tree
// In this method we create a temporary array count[] of size equal to the height of tree. We initialize all values in count as 0. We traverse the tree using preorder traversal and fill the entries in count so that the count array contains count of nodes at each level in Binary Tree.

int getMaxWidth(struct node* root)
{
    int width;
    int h = height(root);

    // Create an array that will store count of nodes at each level
    int *count = (int *)calloc(sizeof(int), h);

     int level = 0;

     // Fill the count array using preorder traversal
     getMaxWidthRecur(root, count, level);

     // Return the maximum value from count array
     return getMax(count, h);
}

// A function that fills count array with count of nodes at every
// level of given binary tree
void getMaxWidthRecur(struct node *root, int count[], int level)
{
     if(root)
    {
         count[level]++;
         getMaxWidthRecur(root->left, count, level+1);
        getMaxWidthRecur(root->right, count, level+1);
    }
}
// UTILITY FUNCTIONS 
// Compute the "height" of a tree -- the number of
//nodes along the longest path from the root node
//down to the farthest leaf node.
//int height(struct node* node)
{
     if (node==NULL)
         return 0;
    else
    {
         /* compute the height of each subtree */
         int lHeight = height(node->left);
         int rHeight = height(node->right);
         /* use the larger one */

        return (lHeight > rHeight)? (lHeight+1): (rHeight+1);
    }
}











//Diameter of a binary tree

int diameter(struct node * tree)
{
    /* base case where tree is empty */
    if (tree == 0)
         return 0; 
     /* get the height of left and right sub-trees */
     int lheight = maxDepth(tree->left);
     int rheight = maxDepth(tree->right);

      /* get the diameter of left and right sub-trees */
      int ldiameter = diameter(tree->left);
      int rdiameter = diameter(tree->right);

      /* Return max of following three
      1) Diameter of left subtree
      2) Diameter of right subtree
      3) Height of left subtree + height of right subtree + 1 */
      return max(lheight + rheight + 1, max(ldiameter, rdiameter));
}
// Time Complexity:O(n^2)



// Foldable Binary Tree
// A tree can be folded if left and right subtrees of the tree are structure wise mirror image of each other. An empty tree is considered as foldable.

bool IsFoldable(struct node *root)
{
     if (root == NULL)
         {  return true;  }

     return IsFoldableUtil(root->left, root->right);
}

/* A utility function that checks if trees with roots as n1 and n2
 are mirror of each other */
bool IsFoldableUtil(struct node *n1, struct node *n2)
{
    /* If both left and right subtrees are NULL,
      then return true */
     if (n1 == NULL && n2 == NULL)
          {  return true;  }

     /* If one of the trees is NULL and other is not,
      then return false */
      if (n1 == NULL || n2 == NULL)
             {  return false; }

      /* Otherwise check if left and right subtrees are mirrors of
       their counterparts */
    return IsFoldableUtil(n1->left, n2->right) && IsFoldableUtil(n1->right, n2->left);
}




// C++ Program to find diagonal 
// sum in a Binary Tree 
#include <iostream> 
#include <stdlib.h> 
#include <map> 
using namespace std; 
  
struct Node 
{ 
    int data; 
    struct Node* left; 
    struct Node* right; 
}; 
  
struct Node* newNode(int data) 
{ 
    struct Node* Node = 
            (struct Node*)malloc(sizeof(struct Node)); 
      
    Node->data = data; 
    Node->left = Node->right = NULL; 
  
    return Node; 
} 
  
// root - root of the binary tree 
// vd - vertical distance diagonally 
// diagonalSum - map to store Diagonal  
// Sum(Passed by Reference) 
void diagonalSumUtil(struct Node* root, 
                int vd, map<int, int> &diagonalSum) 
{ 
    if(!root) 
        return; 
          
    diagonalSum[vd] += root->data; 
  
    // increase the vertical distance if left child 
    diagonalSumUtil(root->left, vd + 1, diagonalSum); 
  
    // vertical distance remains same for right child 
    diagonalSumUtil(root->right, vd, diagonalSum); 
} 
  
// Function to calculate diagonal  
// sum of given binary tree 
void diagonalSum(struct Node* root) 
{ 
  
    // create a map to store Diagonal Sum 
    map<int, int> diagonalSum;  
      
    diagonalSumUtil(root, 0, diagonalSum); 
  
    map<int, int>::iterator it; 
        cout << "Diagonal sum in a binary tree is - "; 
      
    for(it = diagonalSum.begin(); 
                it != diagonalSum.end(); ++it) 
    { 
        cout << it->second << " "; 
    } 
} 
  
// Driver code 
int main() 
{ 
    struct Node* root = newNode(1); 
    root->left = newNode(2); 
    root->right = newNode(3); 
    root->left->left = newNode(9); 
    root->left->right = newNode(6); 
    root->right->left = newNode(4); 
    root->right->right = newNode(5); 
    root->right->left->right = newNode(7); 
    root->right->left->left = newNode(12); 
    root->left->right->left = newNode(11); 
    root->left->left->right = newNode(10); 
  
    diagonalSum(root); 
  
    return 0; 
} 


// A C++ program to find bridges in a given undirected graph 
#include<iostream> 
#include <list> 
#define NIL -1 
using namespace std; 
  
// A class that represents an undirected graph 
class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // A dynamic array of adjacency lists 
    void bridgeUtil(int v, bool visited[], int disc[], int low[], 
                    int parent[]); 
public: 
    Graph(int V);   // Constructor 
    void addEdge(int v, int w);   // to add an edge to graph 
    void bridge();    // prints all bridges 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); 
    adj[w].push_back(v);  // Note: the graph is undirected 
} 
  
// A recursive function that finds and prints bridges using 
// DFS traversal 
// u --> The vertex to be visited next 
// visited[] --> keeps tract of visited vertices 
// disc[] --> Stores discovery times of visited vertices 
// parent[] --> Stores parent vertices in DFS tree 
void Graph::bridgeUtil(int u, bool visited[], int disc[],  
                                  int low[], int parent[]) 
{ 
    // A static variable is used for simplicity, we can  
    // avoid use of static variable by passing a pointer. 
    static int time = 0; 
  
    // Mark the current node as visited 
    visited[u] = true; 
  
    // Initialize discovery time and low value 
    disc[u] = low[u] = ++time; 
  
    // Go through all vertices aadjacent to this 
    list<int>::iterator i; 
    for (i = adj[u].begin(); i != adj[u].end(); ++i) 
    { 
        int v = *i;  // v is current adjacent of u 
  
        // If v is not visited yet, then recur for it 
        if (!visited[v]) 
        { 
            parent[v] = u; 
            bridgeUtil(v, visited, disc, low, parent); 
  
            // Check if the subtree rooted with v has a  
            // connection to one of the ancestors of u 
            low[u]  = min(low[u], low[v]); 
  
            // If the lowest vertex reachable from subtree  
            // under v is  below u in DFS tree, then u-v  
            // is a bridge 
            if (low[v] > disc[u]) 
              cout << u <<" " << v << endl; 
        } 
  
        // Update low value of u for parent function calls. 
        else if (v != parent[u]) 
            low[u]  = min(low[u], disc[v]); 
    } 
} 
  
// DFS based function to find all bridges. It uses recursive  
// function bridgeUtil() 
void Graph::bridge() 
{ 
    // Mark all the vertices as not visited 
    bool *visited = new bool[V]; 
    int *disc = new int[V]; 
    int *low = new int[V]; 
    int *parent = new int[V]; 
  
    // Initialize parent and visited arrays 
    for (int i = 0; i < V; i++) 
    { 
        parent[i] = NIL; 
        visited[i] = false; 
    } 
  
    // Call the recursive helper function to find Bridges 
    // in DFS tree rooted with vertex 'i' 
    for (int i = 0; i < V; i++) 
        if (visited[i] == false) 
            bridgeUtil(i, visited, disc, low, parent); 
} 
  
// Driver program to test above function 
int main() 
{ 
    // Create graphs given in above diagrams 
    cout << "\nBridges in first graph \n"; 
    Graph g1(5); 
    g1.addEdge(1, 0); 
    g1.addEdge(0, 2); 
    g1.addEdge(2, 1); 
    g1.addEdge(0, 3); 
    g1.addEdge(3, 4); 
    g1.bridge(); 
  
    cout << "\nBridges in second graph \n"; 
    Graph g2(4); 
    g2.addEdge(0, 1); 
    g2.addEdge(1, 2); 
    g2.addEdge(2, 3); 
    g2.bridge(); 
  
    cout << "\nBridges in third graph \n"; 
    Graph g3(7); 
    g3.addEdge(0, 1); 
    g3.addEdge(1, 2); 
    g3.addEdge(2, 0); 
    g3.addEdge(1, 3); 
    g3.addEdge(1, 4); 
    g3.addEdge(1, 6); 
    g3.addEdge(3, 5); 
    g3.addEdge(4, 5); 
    g3.bridge(); 
  
    return 0; 
} 



// / CPP program to count subarrays having product 
// less than k. 
#include <iostream> 
#include <vector> 
using namespace std; 
   
int countSubArrayProductLessThanK(const vector<int>& a,  
                                           long long k) 
{ 
    const int n = a.size();     
    long long p = 1; 
    int res = 0; 
    for (int start = 0, end = 0; end < n; end++) { 
  
        // Move right bound by 1 step. Update the product. 
        p *= a[end]; 
          
        // Move left bound so guarantee that p is again  
        // less than k. 
        while (start < end && p >= k)  
            p /= a[start++];         
          
        // If p is less than k, update the counter. 
        // Note that this is working even for (start == end): 
        // it means that the previous window cannot grow  
        // anymore and a single array element is the only  
        // addendum. 
        if (p < k) { 
            int len = end-start+1; 
            res += len; 
        } 
    } 
  
    return res; 
} 
   
// Driver Function to count number of 
// such arrays 
int main() 
{ 
    cout << countSubArrayProductLessThanK({1, 2, 3, 4}, 10) 
         << endl; 
    cout << countSubArrayProductLessThanK({1, 9, 2, 8, 6,  
            4, 3}, 100) << endl; 
    cout << countSubArrayProductLessThanK({5, 3, 2}, 16)  
         << endl; 
    cout << countSubArrayProductLessThanK({100, 200}, 100)  
         << endl; 
    cout << countSubArrayProductLessThanK({100, 200}, 101) 
         << endl; 
} 


// A C++ program to print topological sorting of a graph 
// using indegrees. 
#include<bits/stdc++.h> 
using namespace std; 
  
// Class to represent a graph 
class Graph 
{ 
    int V;    // No. of vertices' 
  
    // Pointer to an array containing adjacency listsList 
    list<int> *adj; 
  
public: 
    Graph(int V);   // Constructor 
  
    // function to add an edge to graph 
    void addEdge(int u, int v); 
  
    // prints a Topological Sort of the complete graph 
    void topologicalSort(); 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int u, int v) 
{ 
    adj[u].push_back(v); 
} 
  
// The function to do Topological Sort. 
void Graph::topologicalSort() 
{ 
    // Create a vector to store indegrees of all 
    // vertices. Initialize all indegrees as 0. 
    vector<int> in_degree(V, 0); 
  
    // Traverse adjacency lists to fill indegrees of 
    // vertices.  This step takes O(V+E) time 
    for (int u=0; u<V; u++) 
    { 
        list<int>::iterator itr; 
        for (itr = adj[u].begin(); itr != adj[u].end(); itr++) 
             in_degree[*itr]++; 
    } 
  
    // Create an queue and enqueue all vertices with 
    // indegree 0 
    queue<int> q; 
    for (int i = 0; i < V; i++) 
        if (in_degree[i] == 0) 
            q.push(i); 
  
    // Initialize count of visited vertices 
    int cnt = 0; 
  
    // Create a vector to store result (A topological 
    // ordering of the vertices) 
    vector <int> top_order; 
  
    // One by one dequeue vertices from queue and enqueue 
    // adjacents if indegree of adjacent becomes 0 
    while (!q.empty()) 
    { 
        // Extract front of queue (or perform dequeue) 
        // and add it to topological order 
        int u = q.front(); 
        q.pop(); 
        top_order.push_back(u); 
  
        // Iterate through all its neighbouring nodes 
        // of dequeued node u and decrease their in-degree 
        // by 1 
        list<int>::iterator itr; 
        for (itr = adj[u].begin(); itr != adj[u].end(); itr++) 
  
            // If in-degree becomes zero, add it to queue 
            if (--in_degree[*itr] == 0) 
                q.push(*itr); 
  
        cnt++; 
    } 
  
    // Check if there was a cycle 
    if (cnt != V) 
    { 
        cout << "There exists a cycle in the graph\n"; 
        return; 
    } 
  
    // Print topological order 
    for (int i=0; i<top_order.size(); i++) 
        cout << top_order[i] << " "; 
    cout << endl; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Create a graph given in the above diagram 
    Graph g(6); 
    g.addEdge(5, 2); 
    g.addEdge(5, 0); 
    g.addEdge(4, 0); 
    g.addEdge(4, 1); 
    g.addEdge(2, 3); 
    g.addEdge(3, 1); 
  
    cout << "Following is a Topological Sort of\n"; 
    g.topologicalSort(); 
  
    return 0; 
} 



// A C++ program to print topological sorting of a DAG 
#include<iostream> 
#include <list> 
#include <stack> 
using namespace std; 
  
// Class to represent a graph 
class Graph 
{ 
    int V;    // No. of vertices' 
  
    // Pointer to an array containing adjacency listsList 
    list<int> *adj; 
  
    // A function used by topologicalSort 
    void topologicalSortUtil(int v, bool visited[], stack<int> &Stack); 
public: 
    Graph(int V);   // Constructor 
  
     // function to add an edge to graph 
    void addEdge(int v, int w); 
  
    // prints a Topological Sort of the complete graph 
    void topologicalSort(); 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); // Add w to v’s list. 
} 
  
// A recursive function used by topologicalSort 
void Graph::topologicalSortUtil(int v, bool visited[],  
                                stack<int> &Stack) 
{ 
    // Mark the current node as visited. 
    visited[v] = true; 
  
    // Recur for all the vertices adjacent to this vertex 
    list<int>::iterator i; 
    for (i = adj[v].begin(); i != adj[v].end(); ++i) 
        if (!visited[*i]) 
            topologicalSortUtil(*i, visited, Stack); 
  
    // Push current vertex to stack which stores result 
    Stack.push(v); 
} 
  
// The function to do Topological Sort. It uses recursive  
// topologicalSortUtil() 
void Graph::topologicalSort() 
{ 
    stack<int> Stack; 
  
    // Mark all the vertices as not visited 
    bool *visited = new bool[V]; 
    for (int i = 0; i < V; i++) 
        visited[i] = false; 
  
    // Call the recursive helper function to store Topological 
    // Sort starting from all vertices one by one 
    for (int i = 0; i < V; i++) 
      if (visited[i] == false) 
        topologicalSortUtil(i, visited, Stack); 
  
    // Print contents of stack 
    while (Stack.empty() == false) 
    { 
        cout << Stack.top() << " "; 
        Stack.pop(); 
    } 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Create a graph given in the above diagram 
    Graph g(6); 
    g.addEdge(5, 2); 
    g.addEdge(5, 0); 
    g.addEdge(4, 0); 
    g.addEdge(4, 1); 
    g.addEdge(2, 3); 
    g.addEdge(3, 1); 
  
    cout << "Following is a Topological Sort of the given graph \n"; 
    g.topologicalSort(); 
  
    return 0; 
} 


Time Complexity: The above algorithm is simply DFS with an extra stack. So time complexity is the same as DFS which is O(V+E).

Note : Here, we can also use vector instead of stack. If the vector is used then print the elements in reverse order to get the topological sorting.

Applications:
Topological Sorting is mainly used for scheduling jobs from the given dependencies among jobs. In computer science, applications of this type arise in instruction scheduling, ordering of formula cell evaluation when recomputing formula values in spreadsheets, logic synthesis, determining the order of compilation tasks to perform in makefiles, data serialization, and resolving symbol dependencies in linkers [2].








// A C++ program to find bridges in a given undirected graph 
#include<iostream> 
#include <list> 
#define NIL -1 
using namespace std; 
  
// A class that represents an undirected graph 
class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // A dynamic array of adjacency lists 
    void bridgeUtil(int v, bool visited[], int disc[], int low[], 
                    int parent[]); 
public: 
    Graph(int V);   // Constructor 
    void addEdge(int v, int w);   // to add an edge to graph 
    void bridge();    // prints all bridges 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); 
    adj[w].push_back(v);  // Note: the graph is undirected 
} 
  
// A recursive function that finds and prints bridges using 
// DFS traversal 
// u --> The vertex to be visited next 
// visited[] --> keeps tract of visited vertices 
// disc[] --> Stores discovery times of visited vertices 
// parent[] --> Stores parent vertices in DFS tree 
void Graph::bridgeUtil(int u, bool visited[], int disc[],  
                                  int low[], int parent[]) 
{ 
    // A static variable is used for simplicity, we can  
    // avoid use of static variable by passing a pointer. 
    static int time = 0; 
  
    // Mark the current node as visited 
    visited[u] = true; 
  
    // Initialize discovery time and low value 
    disc[u] = low[u] = ++time; 
  
    // Go through all vertices aadjacent to this 
    list<int>::iterator i; 
    for (i = adj[u].begin(); i != adj[u].end(); ++i) 
    { 
        int v = *i;  // v is current adjacent of u 
  
        // If v is not visited yet, then recur for it 
        if (!visited[v]) 
        { 
            parent[v] = u; 
            bridgeUtil(v, visited, disc, low, parent); 
  
            // Check if the subtree rooted with v has a  
            // connection to one of the ancestors of u 
            low[u]  = min(low[u], low[v]); 
  
            // If the lowest vertex reachable from subtree  
            // under v is  below u in DFS tree, then u-v  
            // is a bridge 
            if (low[v] > disc[u]) 
              cout << u <<" " << v << endl; 
        } 
  
        // Update low value of u for parent function calls. 
        else if (v != parent[u]) 
            low[u]  = min(low[u], disc[v]); 
    } 
} 
  
// DFS based function to find all bridges. It uses recursive  
// function bridgeUtil() 
void Graph::bridge() 
{ 
    // Mark all the vertices as not visited 
    bool *visited = new bool[V]; 
    int *disc = new int[V]; 
    int *low = new int[V]; 
    int *parent = new int[V]; 
  
    // Initialize parent and visited arrays 
    for (int i = 0; i < V; i++) 
    { 
        parent[i] = NIL; 
        visited[i] = false; 
    } 
  
    // Call the recursive helper function to find Bridges 
    // in DFS tree rooted with vertex 'i' 
    for (int i = 0; i < V; i++) 
        if (visited[i] == false) 
            bridgeUtil(i, visited, disc, low, parent); 
} 
  
// Driver program to test above function 
int main() 
{ 
    // Create graphs given in above diagrams 
    cout << "\nBridges in first graph \n"; 
    Graph g1(5); 
    g1.addEdge(1, 0); 
    g1.addEdge(0, 2); 
    g1.addEdge(2, 1); 
    g1.addEdge(0, 3); 
    g1.addEdge(3, 4); 
    g1.bridge(); 
  
    cout << "\nBridges in second graph \n"; 
    Graph g2(4); 
    g2.addEdge(0, 1); 
    g2.addEdge(1, 2); 
    g2.addEdge(2, 3); 
    g2.bridge(); 
  
    cout << "\nBridges in third graph \n"; 
    Graph g3(7); 
    g3.addEdge(0, 1); 
    g3.addEdge(1, 2); 
    g3.addEdge(2, 0); 
    g3.addEdge(1, 3); 
    g3.addEdge(1, 4); 
    g3.addEdge(1, 6); 
    g3.addEdge(3, 5); 
    g3.addEdge(4, 5); 
    g3.bridge(); 
  
    return 0; 
} 




Cut vertices

// A C++ program to find articulation points in an undirected graph 
#include<iostream> 
#include <list> 
#define NIL -1 
using namespace std; 
  
// A class that represents an undirected graph 
class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // A dynamic array of adjacency lists 
    void APUtil(int v, bool visited[], int disc[], int low[],  
                int parent[], bool ap[]); 
public: 
    Graph(int V);   // Constructor 
    void addEdge(int v, int w);   // function to add an edge to graph 
    void AP();    // prints articulation points 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); 
    adj[w].push_back(v);  // Note: the graph is undirected 
} 
  
// A recursive function that find articulation points using DFS traversal 
// u --> The vertex to be visited next 
// visited[] --> keeps tract of visited vertices 
// disc[] --> Stores discovery times of visited vertices 
// parent[] --> Stores parent vertices in DFS tree 
// ap[] --> Store articulation points 
void Graph::APUtil(int u, bool visited[], int disc[],  
                                      int low[], int parent[], bool ap[]) 
{ 
    // A static variable is used for simplicity, we can avoid use of static 
    // variable by passing a pointer. 
    static int time = 0; 
  
    // Count of children in DFS Tree 
    int children = 0; 
  
    // Mark the current node as visited 
    visited[u] = true; 
  
    // Initialize discovery time and low value 
    disc[u] = low[u] = ++time; 
  
    // Go through all vertices aadjacent to this 
    list<int>::iterator i; 
    for (i = adj[u].begin(); i != adj[u].end(); ++i) 
    { 
        int v = *i;  // v is current adjacent of u 
  
        // If v is not visited yet, then make it a child of u 
        // in DFS tree and recur for it 
        if (!visited[v]) 
        { 
            children++; 
            parent[v] = u; 
            APUtil(v, visited, disc, low, parent, ap); 
  
            // Check if the subtree rooted with v has a connection to 
            // one of the ancestors of u 
            low[u]  = min(low[u], low[v]); 
  
            // u is an articulation point in following cases 
  
            // (1) u is root of DFS tree and has two or more chilren. 
            if (parent[u] == NIL && children > 1) 
               ap[u] = true; 
  
            // (2) If u is not root and low value of one of its child is more 
            // than discovery value of u. 
            if (parent[u] != NIL && low[v] >= disc[u]) 
               ap[u] = true; 
        } 
  
        // Update low value of u for parent function calls. 
        else if (v != parent[u]) 
            low[u]  = min(low[u], disc[v]); 
    } 
} 
  
// The function to do DFS traversal. It uses recursive function APUtil() 
void Graph::AP() 
{ 
    // Mark all the vertices as not visited 
    bool *visited = new bool[V]; 
    int *disc = new int[V]; 
    int *low = new int[V]; 
    int *parent = new int[V]; 
    bool *ap = new bool[V]; // To store articulation points 
  
    // Initialize parent and visited, and ap(articulation point) arrays 
    for (int i = 0; i < V; i++) 
    { 
        parent[i] = NIL; 
        visited[i] = false; 
        ap[i] = false; 
    } 
  
    // Call the recursive helper function to find articulation points 
    // in DFS tree rooted with vertex 'i' 
    for (int i = 0; i < V; i++) 
        if (visited[i] == false) 
            APUtil(i, visited, disc, low, parent, ap); 
  
    // Now ap[] contains articulation points, print them 
    for (int i = 0; i < V; i++) 
        if (ap[i] == true) 
            cout << i << " "; 
} 
  
// Driver program to test above function 
int main() 
{ 
    // Create graphs given in above diagrams 
    cout << "\nArticulation points in first graph \n"; 
    Graph g1(5); 
    g1.addEdge(1, 0); 
    g1.addEdge(0, 2); 
    g1.addEdge(2, 1); 
    g1.addEdge(0, 3); 
    g1.addEdge(3, 4); 
    g1.AP(); 
  
    cout << "\nArticulation points in second graph \n"; 
    Graph g2(4); 
    g2.addEdge(0, 1); 
    g2.addEdge(1, 2); 
    g2.addEdge(2, 3); 
    g2.AP(); 
  
    cout << "\nArticulation points in third graph \n"; 
    Graph g3(7); 
    g3.addEdge(0, 1); 
    g3.addEdge(1, 2); 
    g3.addEdge(2, 0); 
    g3.addEdge(1, 3); 
    g3.addEdge(1, 4); 
    g3.addEdge(1, 6); 
    g3.addEdge(3, 5); 
    g3.addEdge(4, 5); 
    g3.AP(); 
  
    return 0; 
} 

Time Complexity: The above function is simple DFS with additional arrays. 
So time complexity is same as DFS which is O(V+E) for adjacency list representation of graph.




LCA

// C++ Program for Lowest Common Ancestor in a Binary Tree 
// A O(n) solution to find LCA of two given values n1 and n2 
#include <iostream> 
#include <vector> 
  
using namespace std; 
  
// A Binary Tree node 
struct Node 
{ 
    int key; 
    struct Node *left, *right; 
}; 
  
// Utility function creates a new binary tree node with given key 
Node * newNode(int k) 
{ 
    Node *temp = new Node; 
    temp->key = k; 
    temp->left = temp->right = NULL; 
    return temp; 
} 
  
// Finds the path from root node to given root of the tree, Stores the 
// path in a vector path[], returns true if path exists otherwise false 
bool findPath(Node *root, vector<int> &path, int k) 
{ 
    // base case 
    if (root == NULL) return false; 
  
    // Store this node in path vector. The node will be removed if 
    // not in path from root to k 
    path.push_back(root->key); 
  
    // See if the k is same as root's key 
    if (root->key == k) 
        return true; 
  
    // Check if k is found in left or right sub-tree 
    if ( (root->left && findPath(root->left, path, k)) || 
         (root->right && findPath(root->right, path, k)) ) 
        return true; 
  
    // If not present in subtree rooted with root, remove root from 
    // path[] and return false 
    path.pop_back(); 
    return false; 
} 
  
// Returns LCA if node n1, n2 are present in the given binary tree, 
// otherwise return -1 
int findLCA(Node *root, int n1, int n2) 
{ 
    // to store paths to n1 and n2 from the root 
    vector<int> path1, path2; 
  
    // Find paths from root to n1 and root to n1. If either n1 or n2 
    // is not present, return -1 
    if ( !findPath(root, path1, n1) || !findPath(root, path2, n2)) 
          return -1; 
  
    /* Compare the paths to get the first different value */
    int i; 
    for (i = 0; i < path1.size() && i < path2.size() ; i++) 
        if (path1[i] != path2[i]) 
            break; 
    return path1[i-1]; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Let us create the Binary Tree shown in above diagram. 
    Node * root = newNode(1); 
    root->left = newNode(2); 
    root->right = newNode(3); 
    root->left->left = newNode(4); 
    root->left->right = newNode(5); 
    root->right->left = newNode(6); 
    root->right->right = newNode(7); 
    cout << "LCA(4, 5) = " << findLCA(root, 4, 5); 
    cout << "nLCA(4, 6) = " << findLCA(root, 4, 6); 
    cout << "nLCA(3, 4) = " << findLCA(root, 3, 4); 
    cout << "nLCA(2, 4) = " << findLCA(root, 2, 4); 
    return 0; 
} 




Method 2 (Using Single Traversal)
The method 1 finds LCA in O(n) time, but requires three tree traversals plus extra spaces for path arrays. 
If we assume that the keys n1 and n2 are present in Binary Tree, we can find LCA using single traversal of Binary Tree and 
without extra storage for path arrays.
The idea is to traverse the tree starting from root. If any of the given keys (n1 and n2) matches with root,
then root is LCA (assuming that both keys are present). If root doesn’t match with any of the keys, 
we recur for left and right subtree. 
The node which has one key present in its left subtree and the other key present in right subtree is the LCA. 
If both keys lie in left subtree, then left subtree has LCA also, otherwise LCA lies in right subtree.

/* C++ Program to find LCA of n1 and n2 using one traversal of Binary Tree */
#include <iostream> 
  
using namespace std; 
  
// A Binary Tree Node 
struct Node 
{ 
    struct Node *left, *right; 
    int key; 
}; 
  
// Utility function to create a new tree Node 
Node* newNode(int key) 
{ 
    Node *temp = new Node; 
    temp->key = key; 
    temp->left = temp->right = NULL; 
    return temp; 
} 
  
// This function returns pointer to LCA of two given values n1 and n2. 
// This function assumes that n1 and n2 are present in Binary Tree 
struct Node *findLCA(struct Node* root, int n1, int n2) 
{ 
    // Base case 
    if (root == NULL) return NULL; 
  
    // If either n1 or n2 matches with root's key, report 
    // the presence by returning root (Note that if a key is 
    // ancestor of other, then the ancestor key becomes LCA 
    if (root->key == n1 || root->key == n2) 
        return root; 
  
    // Look for keys in left and right subtrees 
    Node *left_lca  = findLCA(root->left, n1, n2); 
    Node *right_lca = findLCA(root->right, n1, n2); 
  
    // If both of the above calls return Non-NULL, then one key 
    // is present in once subtree and other is present in other, 
    // So this node is the LCA 
    if (left_lca && right_lca)  return root; 
  
    // Otherwise check if left subtree or right subtree is LCA 
    return (left_lca != NULL)? left_lca: right_lca; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Let us create binary tree given in the above example 
    Node * root = newNode(1); 
    root->left = newNode(2); 
    root->right = newNode(3); 
    root->left->left = newNode(4); 
    root->left->right = newNode(5); 
    root->right->left = newNode(6); 
    root->right->right = newNode(7); 
    cout << "LCA(4, 5) = " << findLCA(root, 4, 5)->key; 
    cout << "nLCA(4, 6) = " << findLCA(root, 4, 6)->key; 
    cout << "nLCA(3, 4) = " << findLCA(root, 3, 4)->key; 
    cout << "nLCA(2, 4) = " << findLCA(root, 2, 4)->key; 
    return 0; 
}











// C++ program to find lowest common ancestor using parent pointer 
#include <bits/stdc++.h> 
using namespace std; 
  
// A Tree Node 
struct Node 
{ 
    Node *left, *right, *parent; 
    int key; 
}; 
  
// A utility function to create a new BST node 
Node *newNode(int item) 
{ 
    Node *temp = new Node; 
    temp->key = item; 
    temp->parent = temp->left = temp->right = NULL; 
    return temp; 
} 
  
/* A utility function to insert a new node with 
   given key in Binary Search Tree */
Node *insert(Node *node, int key) 
{ 
    /* If the tree is empty, return a new node */
    if (node == NULL) return newNode(key); 
  
    /* Otherwise, recur down the tree */
    if (key < node->key) 
    { 
        node->left  = insert(node->left, key); 
        node->left->parent = node; 
    } 
    else if (key > node->key) 
    { 
        node->right = insert(node->right, key); 
        node->right->parent = node; 
    } 
  
    /* return the (unchanged) node pointer */
    return node; 
} 
  
// To find LCA of nodes n1 and n2 in Binary Tree 
Node *LCA(Node *n1, Node *n2) 
{ 
   // Creata a map to store ancestors of n1 
   map <Node *, bool> ancestors; 
  
   // Insert n1 and all its ancestors in map 
   while (n1 != NULL) 
   { 
       ancestors[n1] = true; 
       n1 = n1->parent; 
   } 
  
   // Check if n2 or any of its ancestors is in 
   // map. 
   while (n2 != NULL) 
   { 
       if (ancestors.find(n2) != ancestors.end()) 
           return n2; 
       n2 = n2->parent; 
   } 
  
   return NULL; 
} 
  
// Driver method to test above functions 
int main(void) 
{ 
    Node * root = NULL; 
  
    root = insert(root, 20); 
    root = insert(root, 8); 
    root = insert(root, 22); 
    root = insert(root, 4); 
    root = insert(root, 12); 
    root = insert(root, 10); 
    root = insert(root, 14); 
  
    Node *n1 = root->left->right->left; 
    Node *n2 = root->left; 
    Node *lca = LCA(n1, n2); 
  
    printf("LCA of %d and %d is %d \n", n1->key, n2->key, lca->key); 
  
    return 0; 
} 








// C++ program to find lowest common ancestor using parent pointer 
#include <bits/stdc++.h> 
using namespace std; 
  
// A Tree Node 
struct Node 
{ 
    Node *left, *right, *parent; 
    int key; 
}; 
  
// A utility function to create a new BST node 
Node *newNode(int item) 
{ 
    Node *temp = new Node; 
    temp->key = item; 
    temp->parent = temp->left = temp->right = NULL; 
    return temp; 
} 
  
/* A utility function to insert a new node with 
given key in Binary Search Tree */
Node *insert(Node *node, int key) 
{ 
    /* If the tree is empty, return a new node */
    if (node == NULL) return newNode(key); 
  
    /* Otherwise, recur down the tree */
    if (key < node->key) 
    { 
        node->left = insert(node->left, key); 
        node->left->parent = node; 
    } 
    else if (key > node->key) 
    { 
        node->right = insert(node->right, key); 
        node->right->parent = node; 
    } 
  
    /* return the (unchanged) node pointer */
    return node; 
} 
  
// A utility function to find depth of a node 
// (distance of it from root) 
int depth(Node *node) 
{ 
    int d = -1; 
    while (node) 
    { 
        ++d; 
        node = node->parent; 
    } 
    return d; 
} 
  
// To find LCA of nodes n1 and n2 in Binary Tree 
Node *LCA(Node *n1, Node *n2) 
{ 
    // Find depths of two nodes and differences 
    int d1 = depth(n1), d2 = depth(n2); 
    int diff = d1 - d2; 
  
    // If n2 is deeper, swap n1 and n2 
    if (diff < 0) 
    { 
        Node * temp = n1; 
        n1 = n2; 
        n2 = temp; 
        diff = -diff; 
    } 
  
    // Move n1 up until it reaches the same level as n2 
    while (diff--) 
        n1 = n1->parent; 
  
    // Now n1 and n2 are at same levels 
    while (n1 && n2) 
    { 
        if (n1 == n2) 
            return n1; 
        n1 = n1->parent; 
        n2 = n2->parent; 
    } 
  
    return NULL; 
} 
  
// Driver method to test above functions 
int main(void) 
{ 
    Node * root = NULL; 
  
    root = insert(root, 20); 
    root = insert(root, 8); 
    root = insert(root, 22); 
    root = insert(root, 4); 
    root = insert(root, 12); 
    root = insert(root, 10); 
    root = insert(root, 14); 
  
    Node *n1 = root->left->right->left; 
    Node *n2 = root->right; 
  
    Node *lca = LCA(n1, n2); 
    printf("LCA of %d and %d is %d \n", n1->key, n2->key, lca->key); 
  
    return 0; 
}








/* C++ Program to find LCA of u and v by reducing the problem to RMQ */
#include<bits/stdc++.h> 
#define V 9               // number of nodes in input tree 
  
int euler[2*V - 1];       // For Euler tour sequence 
int level[2*V - 1];       // Level of nodes in tour sequence 
int firstOccurrence[V+1]; // First occurences of nodes in tour 
int ind;                  // Variable to fill-in euler and level arrays 
  
// A Binary Tree node 
struct Node 
{ 
    int key; 
    struct Node *left, *right; 
}; 
  
// Utility function creates a new binary tree node with given key 
Node * newNode(int k) 
{ 
    Node *temp = new Node; 
    temp->key = k; 
    temp->left = temp->right = NULL; 
    return temp; 
} 
  
// log base 2 of x 
int Log2(int x) 
{ 
    int ans = 0 ; 
    while (x>>=1) ans++; 
    return ans ; 
} 
  
/*  A recursive function to get the minimum value in a given range 
     of array indexes. The following are parameters for this function. 
  
    st    --> Pointer to segment tree 
    index --> Index of current node in the segment tree. Initially 
              0 is passed as root is always at index 0 
    ss & se  --> Starting and ending indexes of the segment represented 
                  by current node, i.e., st[index] 
    qs & qe  --> Starting and ending indexes of query range */
int RMQUtil(int index, int ss, int se, int qs, int qe, int *st) 
{ 
    // If segment of this node is a part of given range, then return 
    //  the min of the segment 
    if (qs <= ss && qe >= se) 
        return st[index]; 
  
    // If segment of this node is outside the given range 
    else if (se < qs || ss > qe) 
        return -1; 
  
    // If a part of this segment overlaps with the given range 
    int mid = (ss + se)/2; 
  
    int q1 = RMQUtil(2*index+1, ss, mid, qs, qe, st); 
    int q2 = RMQUtil(2*index+2, mid+1, se, qs, qe, st); 
  
    if (q1==-1) return q2; 
  
    else if (q2==-1) return q1; 
  
    return (level[q1] < level[q2]) ? q1 : q2; 
} 
  
// Return minimum of elements in range from index qs (quey start) to 
// qe (query end).  It mainly uses RMQUtil() 
int RMQ(int *st, int n, int qs, int qe) 
{ 
    // Check for erroneous input values 
    if (qs < 0 || qe > n-1 || qs > qe) 
    { 
        printf("Invalid Input"); 
        return -1; 
    } 
  
    return RMQUtil(0, 0, n-1, qs, qe, st); 
} 
  
// A recursive function that constructs Segment Tree for array[ss..se]. 
// si is index of current node in segment tree st 
void constructSTUtil(int si, int ss, int se, int arr[], int *st) 
{ 
    // If there is one element in array, store it in current node of 
    // segment tree and return 
    if (ss == se)st[si] = ss; 
  
    else
    { 
        // If there are more than one elements, then recur for left and 
        // right subtrees and store the minimum of two values in this node 
        int mid = (ss + se)/2; 
        constructSTUtil(si*2+1, ss, mid, arr, st); 
        constructSTUtil(si*2+2, mid+1, se, arr, st); 
  
        if (arr[st[2*si+1]] < arr[st[2*si+2]]) 
            st[si] = st[2*si+1]; 
        else
            st[si] = st[2*si+2]; 
    } 
} 
  
/* Function to construct segment tree from given array. This function 
   allocates memory for segment tree and calls constructSTUtil() to 
   fill the allocated memory */
int *constructST(int arr[], int n) 
{ 
    // Allocate memory for segment tree 
  
    // Height of segment tree 
    int x = Log2(n)+1; 
  
    // Maximum size of segment tree 
    int max_size = 2*(1<<x) - 1;  //  2*pow(2,x) -1 
  
    int *st = new int[max_size]; 
  
    // Fill the allocated memory st 
    constructSTUtil(0, 0, n-1, arr, st); 
  
    // Return the constructed segment tree 
    return st; 
} 
  
// Recursive version of the Euler tour of T 
void eulerTour(Node *root, int l) 
{ 
    /* if the passed node exists */
    if (root) 
    { 
        euler[ind] = root->key; // insert in euler array 
        level[ind] = l;         // insert l in level array 
        ind++;                  // increment index 
  
        /* if unvisited, mark first occurrence */
        if (firstOccurrence[root->key] == -1) 
            firstOccurrence[root->key] = ind-1; 
  
        /* tour left subtree if exists, and remark euler 
           and level arrays for parent on return */
        if (root->left) 
        { 
            eulerTour(root->left, l+1); 
            euler[ind]=root->key; 
            level[ind] = l; 
            ind++; 
        } 
  
        /* tour right subtree if exists, and remark euler 
           and level arrays for parent on return */
        if (root->right) 
        { 
            eulerTour(root->right, l+1); 
            euler[ind]=root->key; 
            level[ind] = l; 
            ind++; 
        } 
    } 
} 
  
// Returns LCA of nodes n1, n2 (assuming they are 
//  present in the tree) 
int findLCA(Node *root, int u, int v) 
{ 
    /* Mark all nodes unvisited.  Note that the size of 
        firstOccurrence is 1 as node values which vary from 
        1 to 9 are used as indexes */
    memset(firstOccurrence, -1, sizeof(int)*(V+1)); 
  
    /* To start filling euler and level arrays from index 0 */
    ind = 0; 
  
    /* Start Euler tour with root node on level 0 */
    eulerTour(root, 0); 
  
    /* construct segment tree on level array */
    int *st = constructST(level, 2*V-1); 
  
    /* If v before u in Euler tour.  For RMQ to work, first 
       parameter 'u' must be smaller than second 'v' */
    if (firstOccurrence[u]>firstOccurrence[v]) 
       std::swap(u, v); 
  
    // Starting and ending indexes of query range 
    int qs = firstOccurrence[u]; 
    int qe = firstOccurrence[v]; 
  
    // query for index of LCA in tour 
    int index = RMQ(st, 2*V-1, qs, qe); 
  
    /* return LCA node */
    return euler[index]; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Let us create the Binary Tree as shown in the diagram. 
    Node * root = newNode(1); 
    root->left = newNode(2); 
    root->right = newNode(3); 
    root->left->left = newNode(4); 
    root->left->right = newNode(5); 
    root->right->left = newNode(6); 
    root->right->right = newNode(7); 
    root->left->right->left = newNode(8); 
    root->left->right->right = newNode(9); 
  
    int u = 4, v = 9; 
    printf("The LCA of node %d and node %d is node %d.\n",  
            u, v, findLCA(root, u, v)); 
    return 0; 
} 


Note:

We assume that the nodes queried are present in the tree.
We also assumed that if there are V nodes in tree, then keys (or data) of these nodes are in range from 1 to V.
Time complexity:

Euler tour: Number of nodes is V. For a tree, E = V-1. Euler tour (DFS) will take O(V+E) which is O(2*V) which can be written as O(V).
Segment Tree construction : O(n) where n = V + E = 2*V – 1.
Range Minimum query: O(log(n))
Overall this method takes O(n) time for preprocssing, but takes O(Log n) time for query. 
Therefore, it can be useful when we have a single tree on which we want to perform large number of LCA queries 
(Note that LCA is useful for finding shortest path between two nodes of Binary Tree)

Auxiliary Space:

Euler tour array: O(n) where n = 2*V – 1
Node Levels array: O(n)
First Occurrences array: O(V)
Segment Tree: O(n)
Overall: O(n)




// C++ Implementation of Kosaraju's algorithm to print all SCCs 
#include <iostream> 
#include <list> 
#include <stack> 
using namespace std; 
  
class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // An array of adjacency lists 
  
    // Fills Stack with vertices (in increasing order of finishing 
    // times). The top element of stack has the maximum finishing  
    // time 
    void fillOrder(int v, bool visited[], stack<int> &Stack); 
  
    // A recursive function to print DFS starting from v 
    void DFSUtil(int v, bool visited[]); 
public: 
    Graph(int V); 
    void addEdge(int v, int w); 
  
    // The main function that finds and prints strongly connected 
    // components 
    void printSCCs(); 
  
    // Function that returns reverse (or transpose) of this graph 
    Graph getTranspose(); 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
// A recursive function to print DFS starting from v 
void Graph::DFSUtil(int v, bool visited[]) 
{ 
    // Mark the current node as visited and print it 
    visited[v] = true; 
    cout << v << " "; 
  
    // Recur for all the vertices adjacent to this vertex 
    list<int>::iterator i; 
    for (i = adj[v].begin(); i != adj[v].end(); ++i) 
        if (!visited[*i]) 
            DFSUtil(*i, visited); 
} 
  
Graph Graph::getTranspose() 
{ 
    Graph g(V); 
    for (int v = 0; v < V; v++) 
    { 
        // Recur for all the vertices adjacent to this vertex 
        list<int>::iterator i; 
        for(i = adj[v].begin(); i != adj[v].end(); ++i) 
        { 
            g.adj[*i].push_back(v); 
        } 
    } 
    return g; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); // Add w to v’s list. 
} 
  
void Graph::fillOrder(int v, bool visited[], stack<int> &Stack) 
{ 
    // Mark the current node as visited and print it 
    visited[v] = true; 
  
    // Recur for all the vertices adjacent to this vertex 
    list<int>::iterator i; 
    for(i = adj[v].begin(); i != adj[v].end(); ++i) 
        if(!visited[*i]) 
            fillOrder(*i, visited, Stack); 
  
    // All vertices reachable from v are processed by now, push v  
    Stack.push(v); 
} 
  
// The main function that finds and prints all strongly connected  
// components 
void Graph::printSCCs() 
{ 
    stack<int> Stack; 
  
    // Mark all the vertices as not visited (For first DFS) 
    bool *visited = new bool[V]; 
    for(int i = 0; i < V; i++) 
        visited[i] = false; 
  
    // Fill vertices in stack according to their finishing times 
    for(int i = 0; i < V; i++) 
        if(visited[i] == false) 
            fillOrder(i, visited, Stack); 
  
    // Create a reversed graph 
    Graph gr = getTranspose(); 
  
    // Mark all the vertices as not visited (For second DFS) 
    for(int i = 0; i < V; i++) 
        visited[i] = false; 
  
    // Now process all vertices in order defined by Stack 
    while (Stack.empty() == false) 
    { 
        // Pop a vertex from stack 
        int v = Stack.top(); 
        Stack.pop(); 
  
        // Print Strongly connected component of the popped vertex 
        if (visited[v] == false) 
        { 
            gr.DFSUtil(v, visited); 
            cout << endl; 
        } 
    } 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Create a graph given in the above diagram 
    Graph g(5); 
    g.addEdge(1, 0); 
    g.addEdge(0, 2); 
    g.addEdge(2, 1); 
    g.addEdge(0, 3); 
    g.addEdge(3, 4); 
  
    cout << "Following are strongly connected components in "
            "given graph \n"; 
    g.printSCCs(); 
  
    return 0; 
} 





Maximise array sum after taking non-overlapping sub-arrays of length K
Given an integer array arr[] of length N and an integer K, the task is to select some non-overlapping sub-arrays
such that each sub-array is exactly of length K,
no two sub-arrays are adjacent and sum of all the elements of the selected sub-arrays is maximum.



// C++ implementation of the approach 
#include <bits/stdc++.h> 
#define maxLen 10 
using namespace std; 
  
// To store the states of dp 
int dp[maxLen]; 
  
// To check if a given state 
// has been solved 
bool v[maxLen]; 
  
// To store the prefix-sum 
int prefix_sum[maxLen]; 
  
// Function to fill the prefix_sum[] with 
// the prefix sum of the given array 
void findPrefixSum(int arr[], int n) 
{ 
    prefix_sum[0] = arr[0]; 
    for (int i = 1; i < n; i++) 
        prefix_sum[i] = arr[i] + prefix_sum[i - 1]; 
} 
  
// Function to find the maximum sum subsequence 
// such that no two elements are adjacent 
int maxSum(int arr[], int i, int n, int k) 
{ 
    // Base case 
    if (i + k > n) 
        return 0; 
  
    // To check if a state has 
    // been solved 
    if (v[i]) 
        return dp[i]; 
    v[i] = 1; 
  
    int x; 
  
    if (i == 0) 
        x = prefix_sum[k - 1]; 
    else
        x = prefix_sum[i + k - 1] - prefix_sum[i - 1]; 
  
    // Required recurrence relation 
    dp[i] = max(maxSum(arr, i + 1, n, k), 
                x + maxSum(arr, i + k + 1, n, k)); 
  
    // Returning the value 
    return dp[i]; 
} 
  
// Driver code 
int main() 
{ 
    int arr[] = { 1, 3, 7, 6 }; 
    int n = sizeof(arr) / sizeof(int); 
    int k = 1; 
  
    // Finding prefix-sum 
    findPrefixSum(arr, n); 
  
    // Finding the maximum possible sum 
    cout << maxSum(arr, 0, n, k); 
  
    return 0; 
} 










// Program to find if there is a simple path with 
// weight more than k 
#include<bits/stdc++.h> 
using namespace std; 
  
// iPair ==>  Integer Pair 
typedef pair<int, int> iPair; 
  
// This class represents a dipathted graph using 
// adjacency list representation 
class Graph 
{ 
    int V;    // No. of vertices 
  
    // In a weighted graph, we need to store vertex 
    // and weight pair for every edge 
    list< pair<int, int> > *adj; 
    bool pathMoreThanKUtil(int src, int k, vector<bool> &path); 
  
public: 
    Graph(int V);  // Constructor 
  
    // function to add an edge to graph 
    void addEdge(int u, int v, int w); 
    bool pathMoreThanK(int src, int k); 
}; 
  
// Returns true if graph has path more than k length 
bool Graph::pathMoreThanK(int src, int k) 
{ 
    // Create a path array with nothing included 
    // in path 
    vector<bool> path(V, false); 
  
    // Add source vertex to path 
    path[src] = 1; 
  
    return pathMoreThanKUtil(src, k, path); 
} 
  
// Prints shortest paths from src to all other vertices 
bool Graph::pathMoreThanKUtil(int src, int k, vector<bool> &path) 
{ 
    // If k is 0 or negative, return true; 
    if (k <= 0) 
        return true; 
  
    // Get all adjacent vertices of source vertex src and 
    // recursively explore all paths from src. 
    list<iPair>::iterator i; 
    for (i = adj[src].begin(); i != adj[src].end(); ++i) 
    { 
        // Get adjacent vertex and weight of edge 
        int v = (*i).first; 
        int w = (*i).second; 
  
        // If vertex v is already there in path, then 
        // there is a cycle (we ignore this edge) 
        if (path[v] == true) 
            continue; 
  
        // If weight of is more than k, return true 
        if (w >= k) 
            return true; 
  
        // Else add this vertex to path 
        path[v] = true; 
  
        // If this adjacent can provide a path longer 
        // than k, return true. 
        if (pathMoreThanKUtil(v, k-w, path)) 
            return true; 
  
        // Backtrack 
        path[v] = false; 
    } 
  
    // If no adjacent could produce longer path, return 
    // false 
    return false; 
} 
  
// Allocates memory for adjacency list 
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<iPair> [V]; 
} 
  
// Utility function to an edge (u, v) of weight w 
void Graph::addEdge(int u, int v, int w) 
{ 
    adj[u].push_back(make_pair(v, w)); 
    adj[v].push_back(make_pair(u, w)); 
} 
  
// Driver program to test methods of graph class 
int main() 
{ 
    // create the graph given in above fugure 
    int V = 9; 
    Graph g(V); 
  
    //  making above shown graph 
    g.addEdge(0, 1, 4); 
    g.addEdge(0, 7, 8); 
    g.addEdge(1, 2, 8); 
    g.addEdge(1, 7, 11); 
    g.addEdge(2, 3, 7); 
    g.addEdge(2, 8, 2); 
    g.addEdge(2, 5, 4); 
    g.addEdge(3, 4, 9); 
    g.addEdge(3, 5, 14); 
    g.addEdge(4, 5, 10); 
    g.addEdge(5, 6, 2); 
    g.addEdge(6, 7, 1); 
    g.addEdge(6, 8, 6); 
    g.addEdge(7, 8, 7); 
  
    int src = 0; 
    int k = 62; 
    g.pathMoreThanK(src, k)? cout << "Yes\n" : 
                             cout << "No\n"; 
  
    k = 60; 
    g.pathMoreThanK(src, k)? cout << "Yes\n" : 
                             cout << "No\n"; 
  
    return 0; 
} 














// C++ program for implementation of Ford Fulkerson algorithm 
#include <iostream> 
#include <limits.h> 
#include <string.h> 
#include <queue> 
using namespace std; 
  
// Number of vertices in given graph 
#define V 6 
  
/* Returns true if there is a path from source 's' to sink 't' in 
  residual graph. Also fills parent[] to store the path */
bool bfs(int rGraph[V][V], int s, int t, int parent[]) 
{ 
    // Create a visited array and mark all vertices as not visited 
    bool visited[V]; 
    memset(visited, 0, sizeof(visited)); 
  
    // Create a queue, enqueue source vertex and mark source vertex 
    // as visited 
    queue <int> q; 
    q.push(s); 
    visited[s] = true; 
    parent[s] = -1; 
  
    // Standard BFS Loop 
    while (!q.empty()) 
    { 
        int u = q.front(); 
        q.pop(); 
  
        for (int v=0; v<V; v++) 
        { 
            if (visited[v]==false && rGraph[u][v] > 0) 
            { 
                q.push(v); 
                parent[v] = u; 
                visited[v] = true; 
            } 
        } 
    } 
  
    // If we reached sink in BFS starting from source, then return 
    // true, else false 
    return (visited[t] == true); 
} 
  
// Returns the maximum flow from s to t in the given graph 
int fordFulkerson(int graph[V][V], int s, int t) 
{ 
    int u, v; 
  
    // Create a residual graph and fill the residual graph with 
    // given capacities in the original graph as residual capacities 
    // in residual graph 
    int rGraph[V][V]; // Residual graph where rGraph[i][j] indicates  
                     // residual capacity of edge from i to j (if there 
                     // is an edge. If rGraph[i][j] is 0, then there is not)   
    for (u = 0; u < V; u++) 
        for (v = 0; v < V; v++) 
             rGraph[u][v] = graph[u][v]; 
  
    int parent[V];  // This array is filled by BFS and to store path 
  
    int max_flow = 0;  // There is no flow initially 
  
    // Augment the flow while tere is path from source to sink 
    while (bfs(rGraph, s, t, parent)) 
    { 
        // Find minimum residual capacity of the edges along the 
        // path filled by BFS. Or we can say find the maximum flow 
        // through the path found. 
        int path_flow = INT_MAX; 
        for (v=t; v!=s; v=parent[v]) 
        { 
            u = parent[v]; 
            path_flow = min(path_flow, rGraph[u][v]); 
        } 
  
        // update residual capacities of the edges and reverse edges 
        // along the path 
        for (v=t; v != s; v=parent[v]) 
        { 
            u = parent[v]; 
            rGraph[u][v] -= path_flow; 
            rGraph[v][u] += path_flow; 
        } 
  
        // Add path flow to overall flow 
        max_flow += path_flow; 
    } 
  
    // Return the overall flow 
    return max_flow; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    // Let us create a graph shown in the above example 
    int graph[V][V] = { {0, 16, 13, 0, 0, 0}, 
                        {0, 0, 10, 12, 0, 0}, 
                        {0, 4, 0, 0, 14, 0}, 
                        {0, 0, 9, 0, 0, 20}, 
                        {0, 0, 0, 7, 0, 4}, 
                        {0, 0, 0, 0, 0, 0} 
                      }; 
  
    cout << "The maximum possible flow is " << fordFulkerson(graph, 0, 5); 
  
    return 0; 
} 








// C++ program to count Minimum number 
// of jumps to reach end 
#include<bits/stdc++.h> 
using namespace std; 
  
int max(int x, int y) 
 {  
  return (x > y)? x: y;  
 } 
  
// Returns minimum number of jumps to reach arr[n-1] from arr[0] 
int minJumps(int arr[], int n) 
{ 
      
    // The number of jumps needed to reach the starting index is 0 
    if (n <= 1) 
        return 0; 
  
    // Return -1 if not possible to jump 
    if (arr[0] == 0) 
        return -1; 
  
    // initialization 
    int maxReach = arr[0];  // stores all time the maximal reachable index in the array. 
    int step = arr[0];      // stores the number of steps we can still take 
    int jump =1;//stores the number of jumps necessary to reach that maximal reachable position. 
  
    // Start traversing array 
    int i=1; 
    for (i = 1; i < n; i++) 
    { 
        // Check if we have reached the end of the array 
        if (i == n-1) 
            return jump; 
  
        // updating maxReach 
        maxReach = max(maxReach, i+arr[i]); 
  
        // we use a step to get to the current index 
        step--; 
  
        // If no further steps left 
        if (step == 0) 
        { 
            // we must have used a jump 
            jump++; 
  
            // Check if the current index/position or lesser index 
            // is the maximum reach point from the previous indexes 
            if(i >= maxReach) 
                return -1; 
  
            // re-initialize the steps to the amount 
            // of steps to reach maxReach from position i. 
            step = maxReach - i; 
        } 
    } 
  
    return -1; 
} 
  
// Driver program to test above function 
int main() 
{ 
    int arr[]={1, 3, 5, 8, 9, 2, 6, 7, 6, 8, 9}; 
    int size = sizeof(arr)/sizeof(int); 
  
    // Calling the minJumps function 
    cout<<("Minimum number of jumps to reach end is %d ", minJumps(arr,size)); 
    return 0; 
} 
// This code is contributed by 
// Shashank_Sharma





// A C++ Program to detect cycle in an undirected graph 
#include<iostream> 
#include <list> 
#include <limits.h> 
using namespace std; 
  
// Class for an undirected graph 
class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // Pointer to an array containing adjacency lists 
    bool isCyclicUtil(int v, bool visited[], int parent); 
public: 
    Graph(int V);   // Constructor 
    void addEdge(int v, int w);   // to add an edge to graph 
    bool isCyclic();   // returns true if there is a cycle 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); // Add w to v’s list. 
    adj[w].push_back(v); // Add v to w’s list. 
} 
  
// A recursive function that uses visited[] and parent to detect 
// cycle in subgraph reachable from vertex v. 
bool Graph::isCyclicUtil(int v, bool visited[], int parent) 
{ 
    // Mark the current node as visited 
    visited[v] = true; 
  
    // Recur for all the vertices adjacent to this vertex 
    list<int>::iterator i; 
    for (i = adj[v].begin(); i != adj[v].end(); ++i) 
    { 
        // If an adjacent is not visited, then recur for that adjacent 
        if (!visited[*i]) 
        { 
           if (isCyclicUtil(*i, visited, v)) 
              return true; 
        } 
  
        // If an adjacent is visited and not parent of current vertex, 
        // then there is a cycle. 
        else if (*i != parent) 
           return true; 
    } 
    return false; 
} 
  
// Returns true if the graph contains a cycle, else false. 
bool Graph::isCyclic() 
{ 
    // Mark all the vertices as not visited and not part of recursion 
    // stack 
    bool *visited = new bool[V]; 
    for (int i = 0; i < V; i++) 
        visited[i] = false; 
  
    // Call the recursive helper function to detect cycle in different 
    // DFS trees 
    for (int u = 0; u < V; u++) 
        if (!visited[u]) // Don't recur for u if it is already visited 
          if (isCyclicUtil(u, visited, -1)) 
             return true; 
  
    return false; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    Graph g1(5); 
    g1.addEdge(1, 0); 
    g1.addEdge(0, 2); 
    g1.addEdge(2, 1); 
    g1.addEdge(0, 3); 
    g1.addEdge(3, 4); 
    g1.isCyclic()? cout << "Graph contains cycle\n": 
                   cout << "Graph doesn't contain cycle\n"; 
  
    Graph g2(3); 
    g2.addEdge(0, 1); 
    g2.addEdge(1, 2); 
    g2.isCyclic()? cout << "Graph contains cycle\n": 
                   cout << "Graph doesn't contain cycle\n"; 
  
    return 0; 
} 







// C/C++ program to find maximum path sum in Binary Tree 
#include<bits/stdc++.h> 
using namespace std; 
  
// A binary tree node 
struct Node 
{ 
    int data; 
    struct Node* left, *right; 
}; 
  
// A utility function to allocate a new node 
struct Node* newNode(int data) 
{ 
    struct Node* newNode = new Node; 
    newNode->data = data; 
    newNode->left = newNode->right = NULL; 
    return (newNode); 
} 
  
// This function returns overall maximum path sum in 'res' 
// And returns max path sum going through root. 
int findMaxUtil(Node* root, int &res) 
{ 
    //Base Case 
    if (root == NULL) 
        return 0; 
  
    // l and r store maximum path sum going through left and 
    // right child of root respectively 
    int l = findMaxUtil(root->left,res); 
    int r = findMaxUtil(root->right,res); 
  
    // Max path for parent call of root. This path must 
    // include at-most one child of root 
    int max_single = max(max(l, r) + root->data, root->data); 
  
    // Max Top represents the sum when the Node under 
    // consideration is the root of the maxsum path and no 
    // ancestors of root are there in max sum path 
    int max_top = max(max_single, l + r + root->data); 
  
    res = max(res, max_top); // Store the Maximum Result. 
  
    return max_single; 
} 
  
// Returns maximum path sum in tree with given root 
int findMaxSum(Node *root) 
{ 
    // Initialize result 
    int res = INT_MIN; 
  
    // Compute and return result 
    findMaxUtil(root, res); 
    return res; 
} 
  
// Driver program 
int main(void) 
{ 
    struct Node *root = newNode(10); 
    root->left        = newNode(2); 
    root->right       = newNode(10); 
    root->left->left  = newNode(20); 
    root->left->right = newNode(1); 
    root->right->right = newNode(-25); 
    root->right->right->left   = newNode(3); 
    root->right->right->right  = newNode(4); 
    cout << "Max path sum is " << findMaxSum(root); 
    return 0; 
} 





/ We can use stl container list as a double 
// ended queue to store the cache keys, with 
// the descending time of reference from front 
// to back and a set container to check presence 
// of a key. But to fetch the address of the key 
// in the list using find(), it takes O(N) time. 
// This can be optimized by storing a reference 
//     (iterator) to each key in a hash map. 
#include <bits/stdc++.h> 
using namespace std; 
  
class LRUCache { 
    // store keys of cache 
    list<int> dq; 
  
    // store references of key in cache 
    unordered_map<int, list<int>::iterator> ma; 
    int csize; // maximum capacity of cache 
  
public: 
    LRUCache(int); 
    void refer(int); 
    void display(); 
}; 
  
// Declare the size 
LRUCache::LRUCache(int n) 
{ 
    csize = n; 
} 
  
// Refers key x with in the LRU cache 
void LRUCache::refer(int x) 
{ 
    // not present in cache 
    if (ma.find(x) == ma.end()) { 
        // cache is full 
        if (dq.size() == csize) { 
            // delete least recently used element 
            int last = dq.back(); 
  
            // Pops the last elmeent 
            dq.pop_back(); 
  
            // Erase the last 
            ma.erase(last); 
        } 
    } 
  
    // present in cache 
    else
        dq.erase(ma[x]); 
  
    // update reference 
    dq.push_front(x); 
    ma[x] = dq.begin(); 
} 
  
// Function to display contents of cache 
void LRUCache::display() 
{ 
  
    // Iterate in the deque and print 
    // all the elements in it 
    for (auto it = dq.begin(); it != dq.end(); 
         it++) 
        cout << (*it) << " "; 
  
    cout << endl; 
} 
  
// Driver Code 
int main() 
{ 
    LRUCache ca(4); 
  
    ca.refer(1); 
    ca.refer(2); 
    ca.refer(3); 
    ca.refer(1); 
    ca.refer(4); 
    ca.refer(5); 
    ca.display(); 
  
    return 0; 
} 
// This code is contributed by Satish Srinivas 



A Simple Solution is to consider every petrol pumps as a starting point and see if there is a possible tour. If we find a starting point with a feasible solution, we return that starting point. The worst case time complexity of this solution is O(n^2).

An efficient approach is to use a Queue to store the current tour. We first enqueue first petrol pump to the queue, we keep enqueueing petrol pumps till we either complete the tour, or the current amount of petrol becomes negative. If the amount becomes negative, then we keep dequeuing petrol pumps until the queue becomes empty.

Instead of creating a separate queue, we use the given array itself as a queue. We maintain two index variables start and end that represent the rear and front of the queue.

// C++ program to find circular tour for a truck  
#include <bits/stdc++.h> 
using namespace std;  
  
// A petrol pump has petrol and distance to next petrol pump  
class petrolPump  
{ 
    public: 
    int petrol;  
    int distance;  
};  
  
// The function returns starting point if there is a possible solution,  
// otherwise returns -1  
int printTour(petrolPump arr[], int n)  
{  
    // Consider first petrol pump as a starting point  
    int start = 0;  
    int end = 1;  
  
    int curr_petrol = arr[start].petrol - arr[start].distance;  
  
    /* Run a loop while all petrol pumps are not visited.  
    And we have reached first petrol pump again with 0 or more petrol */
    while (end != start || curr_petrol < 0)  
    {  
        // If curremt amount of petrol in truck becomes less than 0, then  
        // remove the starting petrol pump from tour  
        while (curr_petrol < 0 && start != end)  
        {  
            // Remove starting petrol pump. Change start  
            curr_petrol -= arr[start].petrol - arr[start].distance;  
            start = (start + 1) % n;  
  
            // If 0 is being considered as start again, then there is no  
            // possible solution  
            if (start == 0)  
            return -1;  
        }  
  
        // Add a petrol pump to current tour  
        curr_petrol += arr[end].petrol - arr[end].distance;  
  
        end = (end + 1) % n;  
    }  
  
    // Return starting point  
    return start;  
}  
  
// Driver code  
int main()  
{  
    petrolPump arr[] = {{6, 4}, {3, 6}, {7, 3}};  
  
    int n = sizeof(arr)/sizeof(arr[0]);  
    int start = printTour(arr, n);  
  
    (start == -1)? cout<<"No solution": cout<<"Start = "<<start;  
  
    return 0;  
}  
  
  






// Remove all nodes which don’t lie in any path with sum>= k

#include <stdio.h> 
#include <stdlib.h> 
  
// A Binary Tree Node 
struct Node 
{ 
    int data; 
    struct Node *left, *right; 
}; 
  
// A utility function to create a new Binary 
// Tree node with given data 
struct Node* newNode(int data) 
{ 
    struct Node* node = 
       (struct Node*) malloc(sizeof(struct Node)); 
    node->data = data; 
    node->left = node->right = NULL; 
    return node; 
} 
  
// print the tree in LVR (Inorder traversal) way. 
void print(struct Node *root) 
{ 
    if (root != NULL) 
    { 
        print(root->left); 
        printf("%d ",root->data); 
        print(root->right); 
    } 
} 
  
/* Main function which truncates the binary tree. */
struct Node *prune(struct Node *root, int sum) 
{ 
    // Base Case 
    if (root == NULL) return NULL; 
  
    // Recur for left and right subtrees 
    root->left = prune(root->left, sum - root->data); 
    root->right = prune(root->right, sum - root->data); 
  
    // If we reach leaf whose data is smaller than sum, 
    // we delete the leaf.  An important thing to note 
    // is a non-leaf node can become leaf when its 
    // chilren are deleted. 
    if (root->left==NULL && root->right==NULL) 
    { 
        if (root->data < sum) 
        { 
            free(root); 
            return NULL; 
        } 
    } 
  
    return root; 
} 
  
// Driver program to test above function 
int main() 
{ 
    int k = 45; 
    struct Node *root = newNode(1); 
    root->left = newNode(2); 
    root->right = newNode(3); 
    root->left->left = newNode(4); 
    root->left->right = newNode(5); 
    root->right->left = newNode(6); 
    root->right->right = newNode(7); 
    root->left->left->left = newNode(8); 
    root->left->left->right = newNode(9); 
    root->left->right->left = newNode(12); 
    root->right->right->left = newNode(10); 
    root->right->right->left->right = newNode(11); 
    root->left->left->right->left = newNode(13); 
    root->left->left->right->right = newNode(14); 
    root->left->left->right->right->left = newNode(15); 
  
    printf("Tree before truncation\n"); 
    print(root); 
  
    root = prune(root, k); // k is 45 
  
    printf("\n\nTree after truncation\n"); 
    print(root); 
  
    return 0; 
} 







// C++ program to check if binary tree 
// is subtree of another binary tree 
#include<bits/stdc++.h> 
using namespace std; 
  
/* A binary tree node has data,  
left child and right child */
class node  
{  
    public: 
    int data;  
    node* left;  
    node* right;  
};  
  
/* A utility function to check  
whether trees with roots as root1 and  
root2 are identical or not */
bool areIdentical(node * root1, node *root2)  
{  
    /* base cases */
    if (root1 == NULL && root2 == NULL)  
        return true;  
  
    if (root1 == NULL || root2 == NULL)  
        return false;  
  
    /* Check if the data of both roots is  
    same and data of left and right  
    subtrees are also same */
    return (root1->data == root2->data &&  
            areIdentical(root1->left, root2->left) &&  
            areIdentical(root1->right, root2->right) );  
}  
  
  
/* This function returns true if S  
is a subtree of T, otherwise false */
bool isSubtree(node *T, node *S)  
{  
    /* base cases */
    if (S == NULL)  
        return true;  
  
    if (T == NULL)  
        return false;  
  
    /* Check the tree with root as current node */
    if (areIdentical(T, S))  
        return true;  
  
    /* If the tree with root as current  
    node doesn't match then try left  
    and right subtrees one by one */
    return isSubtree(T->left, S) ||  
        isSubtree(T->right, S);  
}  
  
  
/* Helper function that allocates  
a new node with the given data  
and NULL left and right pointers. */
node* newNode(int data)  
{  
    node* Node = new node();  
    Node->data = data;  
    Node->left = NULL;  
    Node->right = NULL;  
    return(Node);  
}  
  
/* Driver code*/
int main()  
{  
    // TREE 1  
    /* Construct the following tree  
            26  
            / \  
        10 3  
        / \ \  
    4 6 3  
    \  
        30  
    */
    node *T = newNode(26);  
    T->right         = newNode(3);  
    T->right->right = newNode(3);  
    T->left         = newNode(10);  
    T->left->left     = newNode(4);  
    T->left->left->right = newNode(30);  
    T->left->right     = newNode(6);  
  
    // TREE 2  
    /* Construct the following tree  
        10  
        / \  
    4 6  
    \  
        30  
    */
    node *S = newNode(10);  
    S->right     = newNode(6);  
    S->left     = newNode(4);  
    S->left->right = newNode(30);  
  
  
    if (isSubtree(T, S))  
        cout << "Tree 2 is subtree of Tree 1";  
    else
        cout << "Tree 2 is not a subtree of Tree 1";  
  
    return 0;  
}  

// check-binary-tree-subtree-another-binary-tree-set-2/
#include <cstring> 
#include <iostream> 
using namespace std; 
#define MAX 100 

// Structure of a tree node 
struct Node { 
	char key; 
	struct Node *left, *right; 
}; 

// A utility function to create a new BST node 
Node* newNode(char item) 
{ 
	Node* temp = new Node; 
	temp->key = item; 
	temp->left = temp->right = NULL; 
	return temp; 
} 

// A utility function to store inorder traversal of tree rooted 
// with root in an array arr[]. Note that i is passed as reference 
void storeInorder(Node* root, char arr[], int& i) 
{ 
	if (root == NULL) { 
		arr[i++] = '$'; 
		return; 
	} 
	storeInorder(root->left, arr, i); 
	arr[i++] = root->key; 
	storeInorder(root->right, arr, i); 
} 

// A utility function to store preorder traversal of tree rooted 
// with root in an array arr[]. Note that i is passed as reference 
void storePreOrder(Node* root, char arr[], int& i) 
{ 
	if (root == NULL) { 
		arr[i++] = '$'; 
		return; 
	} 
	arr[i++] = root->key; 
	storePreOrder(root->left, arr, i); 
	storePreOrder(root->right, arr, i); 
} 

/* This function returns true if S is a subtree of T, otherwise false */
bool isSubtree(Node* T, Node* S) 
{ 
	/* base cases */
	if (S == NULL) 
		return true; 
	if (T == NULL) 
		return false; 

	// Store Inorder traversals of T and S in inT[0..m-1] 
	// and inS[0..n-1] respectively 
	int m = 0, n = 0; 
	char inT[MAX], inS[MAX]; 
	storeInorder(T, inT, m); 
	storeInorder(S, inS, n); 
	inT[m] = '\0', inS[n] = '\0'; 
    cout << inT << endl;

	// If inS[] is not a substring of preS[], return false 
	if (strstr(inT, inS) == NULL) 
		return false; 

	// Store Preorder traversals of T and S in inT[0..m-1] 
	// and inS[0..n-1] respectively 
	m = 0, n = 0; 
	char preT[MAX], preS[MAX]; 
	storePreOrder(T, preT, m); 
	storePreOrder(S, preS, n); 
	preT[m] = '\0', preS[n] = '\0'; 

	// If inS[] is not a substring of preS[], return false 
	// Else return true 
	return (strstr(preT, preS) != NULL); 
} 

// Driver program to test above function 
int main() 
{ 
	Node* T = newNode('a'); 
	T->left = newNode('b'); 
	T->right = newNode('d'); 
	T->left->left = newNode('c'); 
	T->right->right = newNode('e'); 

	Node* S = newNode('a'); 
	S->left = newNode('b'); 
	S->left->left = newNode('c'); 
	S->right = newNode('d'); 

	if (isSubtree(T, S)) 
		cout << "Yes: S is a subtree of T"; 
	else
		cout << "No: S is NOT a subtree of T"; 

	return 0; 
} 

strstr function usage
/* strstr example */
#include <stdio.h>
#include <string.h>

int main ()
{
  char str[] ="This is a simple string";
  char * pch;
  pch = strstr (str,"simple");
  strncpy (pch,"sample",6);
  puts (str);
  return 0;
}



/* C++ implementation to convert infix expression to postfix*/
// Note that here we use std::stack  for Stack operations 
#include<bits/stdc++.h> 
using namespace std; 
  
//Function to return precedence of operators 
int prec(char c) 
{ 
    if(c == '^') 
    return 3; 
    else if(c == '*' || c == '/') 
    return 2; 
    else if(c == '+' || c == '-') 
    return 1; 
    else
    return -1; 
} 
  
// The main function to convert infix expression 
//to postfix expression 
void infixToPostfix(string s) 
{ 
    std::stack<char> st; 
    st.push('N'); 
    int l = s.length(); 
    string ns; 
    for(int i = 0; i < l; i++) 
    { 
        // If the scanned character is an operand, add it to output string. 
        if((s[i] >= 'a' && s[i] <= 'z')||(s[i] >= 'A' && s[i] <= 'Z')) 
        ns+=s[i]; 
  
        // If the scanned character is an ‘(‘, push it to the stack. 
        else if(s[i] == '(') 
          
        st.push('('); 
          
        // If the scanned character is an ‘)’, pop and to output string from the stack 
        // until an ‘(‘ is encountered. 
        else if(s[i] == ')') 
        { 
            while(st.top() != 'N' && st.top() != '(') 
            { 
                char c = st.top(); 
                st.pop(); 
               ns += c; 
            } 
            if(st.top() == '(') 
            { 
                char c = st.top(); 
                st.pop(); 
            } 
        } 
          
        //If an operator is scanned 
        else{ 
            while(st.top() != 'N' && prec(s[i]) <= prec(st.top())) 
            { 
                char c = st.top(); 
                st.pop(); 
                ns += c; 
            } 
            st.push(s[i]); 
        } 
  
    } 
    //Pop all the remaining elements from the stack 
    while(st.top() != 'N') 
    { 
        char c = st.top(); 
        st.pop(); 
        ns += c; 
    } 
      
    cout << ns << endl; 
  
} 
  
//Driver program to test above functions 
int main() 
{ 
    string exp = "a+b*(c^d-e)^(f+g*h)-i"; 
    infixToPostfix(exp); 
    return 0; 
} 





// CPP program to find infix for 
// a given postfix. 
#include <bits/stdc++.h> 
using namespace std; 
  
bool isOperand(char x) 
{ 
   return (x >= 'a' && x <= 'z') || 
          (x >= 'A' && x <= 'Z'); 
} 
  
// Get Infix for a given postfix 
// expression 
string getInfix(string exp) 
{ 
    stack<string> s; 
  
    for (int i=0; exp[i]!='\0'; i++) 
    { 
        // Push operands 
        if (isOperand(exp[i])) 
        { 
           string op(1, exp[i]); 
           s.push(op); 
        } 
  
        // We assume that input is 
        // a valid postfix and expect 
        // an operator. 
        else
        { 
            string op1 = s.top(); 
            s.pop(); 
            string op2 = s.top(); 
            s.pop(); 
            s.push("(" + op2 + exp[i] + 
                   op1 + ")"); 
        } 
    } 
  
    // There must be a single element 
    // in stack now which is the required 
    // infix. 
    return s.top(); 
} 
  
// Driver code 
int main() 
{ 
    string exp = "ab*c+"; 
    cout << getInfix(exp); 
    return 0; 
} 




// CPP program to convert infix to prefix 
#include <bits/stdc++.h> 
using namespace std; 
  
bool isOperator(char c) 
{ 
    return (!isalpha(c) && !isdigit(c)); 
} 
  
int getPriority(char C) 
{ 
    if (C == '-' || C == '+') 
        return 1; 
    else if (C == '*' || C == '/') 
        return 2; 
    else if (C == '^') 
        return 3; 
    return 0; 
} 
  
string infixToPostfix(string infix) 
{ 
    infix = '(' + infix + ')'; 
    int l = infix.size(); 
    stack<char> char_stack; 
    string output; 
  
    for (int i = 0; i < l; i++) { 
  
        // If the scanned character is an  
        // operand, add it to output. 
        if (isalpha(infix[i]) || isdigit(infix[i])) 
            output += infix[i]; 
  
        // If the scanned character is an 
        // ‘(‘, push it to the stack. 
        else if (infix[i] == '(') 
            char_stack.push('('); 
  
        // If the scanned character is an 
        // ‘)’, pop and output from the stack  
        // until an ‘(‘ is encountered. 
        else if (infix[i] == ')') { 
  
            while (char_stack.top() != '(') { 
                output += char_stack.top(); 
                char_stack.pop(); 
            } 
  
            // Remove '(' from the stack 
            char_stack.pop();  
        } 
  
        // Operator found  
        else { 
              
            if (isOperator(char_stack.top())) { 
                while (getPriority(infix[i]) 
                   <= getPriority(char_stack.top())) { 
                    output += char_stack.top(); 
                    char_stack.pop(); 
                } 
  
                // Push current Operator on stack 
                char_stack.push(infix[i]); 
            } 
        } 
    } 
    return output; 
} 
  
string infixToPrefix(string infix) 
{ 
    /* Reverse String 
     * Replace ( with ) and vice versa 
     * Get Postfix 
     * Reverse Postfix  *  */
    int l = infix.size(); 
  
    // Reverse infix 
    reverse(infix.begin(), infix.end()); 
  
    // Replace ( with ) and vice versa 
    for (int i = 0; i < l; i++) { 
  
        if (infix[i] == '(') { 
            infix[i] = ')'; 
            i++; 
        } 
        else if (infix[i] == ')') { 
            infix[i] = '('; 
            i++; 
        } 
    } 
  
    string prefix = infixToPostfix(infix); 
  
    // Reverse postfix 
    reverse(prefix.begin(), prefix.end()); 
  
    return prefix; 
} 
  
// Driver code 
int main() 
{ 
    string s = ("(a-b/c)*(a/k-l)"); 
    cout << infixToPrefix(s) << std::endl; 
    return 0; 
} 





// CPP Program to convert prefix to Infix 
#include <iostream> 
#include <stack> 
using namespace std; 
  
// function to check if character is operator or not 
bool isOperator(char x) { 
  switch (x) { 
  case '+': 
  case '-': 
  case '/': 
  case '*': 
    return true; 
  } 
  return false; 
} 
  
// Convert prefix to Infix expression 
string preToInfix(string pre_exp) { 
  stack<string> s; 
  
  // length of expression 
  int length = pre_exp.size(); 
  
  // reading from right to left 
  for (int i = length - 1; i >= 0; i--) { 
  
    // check if symbol is operator 
    if (isOperator(pre_exp[i])) { 
  
      // pop two operands from stack 
      string op1 = s.top();   s.pop(); 
      string op2 = s.top();   s.pop(); 
  
      // concat the operands and operator 
      string temp = "(" + op1 + pre_exp[i] + op2 + ")"; 
  
      // Push string temp back to stack 
      s.push(temp); 
    } 
  
    // if symbol is an operand 
    else { 
  
      // push the operand to the stack 
      s.push(string(1, pre_exp[i])); 
    } 
  } 
  
  // Stack now contains the Infix expression 
  return s.top(); 
} 
  
// Driver Code 
int main() { 
  string pre_exp = "*-A/BC-/AKL"; 
  cout << "Infix : " << preToInfix(pre_exp); 
  return 0; 
} 






// CPP Program to convert prefix to postfix 
#include <iostream> 
#include <stack> 
using namespace std; 
  
// funtion to check if character is operator or not 
bool isOperator(char x) { 
  switch (x) { 
  case '+': 
  case '-': 
  case '/': 
  case '*': 
    return true; 
  } 
  return false; 
} 
  
// Convert prefix to Postfix expression 
string preToPost(string pre_exp) { 
  
  stack<string> s; 
  
  // length of expression 
  int length = pre_exp.size(); 
  
  // reading from right to left 
  for (int i = length - 1; i >= 0; i--) { 
  
    // check if symbol is operator 
    if (isOperator(pre_exp[i])) { 
  
      // pop two operands from stack 
      string op1 = s.top(); s.pop(); 
      string op2 = s.top(); s.pop(); 
  
      // concat the operands and operator 
      string temp = op1 + op2 + pre_exp[i]; 
  
      // Push string temp back to stack 
      s.push(temp); 
    } 
  
    // if symbol is an operand 
    else { 
  
      // push the operand to the stack 
      s.push(string(1, pre_exp[i])); 
    } 
  } 
  
  // stack contains only the Postfix expression 
  return s.top(); 
} 
  
// Driver Code 
int main() { 
  string pre_exp = "*-A/BC-/AKL"; 
  cout << "Postfix : " << preToPost(pre_exp); 
  return 0; 
} 






// CPP Program to convert postfix to prefix 
#include <iostream> 
#include <stack> 
using namespace std; 
  
// function to check if character is operator or not 
bool isOperator(char x) 
{ 
    switch (x) { 
    case '+': 
    case '-': 
    case '/': 
    case '*': 
        return true; 
    } 
    return false; 
} 
  
// Convert postfix to Prefix expression 
string postToPre(string post_exp) 
{ 
    stack<string> s; 
  
    // length of expression 
    int length = post_exp.size(); 
  
    // reading from right to left 
    for (int i = 0; i < length; i++) { 
  
        // check if symbol is operator 
        if (isOperator(post_exp[i])) { 
  
            // pop two operands from stack 
            string op1 = s.top(); 
            s.pop(); 
            string op2 = s.top(); 
            s.pop(); 
  
            // concat the operands and operator 
            string temp = post_exp[i] + op2 + op1; 
  
            // Push string temp back to stack 
            s.push(temp); 
        } 
  
        // if symbol is an operand 
        else { 
  
            // push the operand to the stack 
            s.push(string(1, post_exp[i])); 
        } 
    } 
  
    // stack[0] contains the Prefix expression 
    return s.top(); 
} 
  
// Driver Code 
int main() 
{ 
    string post_exp = "ABC/-AK/L-*"; 
    cout << "Prefix : " << postToPre(post_exp); 
    return 0; 
} 





// C++ program to find maximum rectangular area in 
// linear time 
#include<iostream> 
#include<stack> 
using namespace std; 
  
// The main function to find the maximum rectangular  
// area under given histogram with n bars 
int getMaxArea(int hist[], int n) 
{ 
    // Create an empty stack. The stack holds indexes  
    // of hist[] array. The bars stored in stack are  
    // always in increasing order of their heights. 
    stack<int> s; 
  
    int max_area = 0; // Initialize max area 
    int tp;  // To store top of stack 
    int area_with_top; // To store area with top bar 
                       // as the smallest bar 
  
    // Run through all bars of given histogram 
    int i = 0; 
    while (i < n) 
    { 
        // If this bar is higher than the bar on top  
        // stack, push it to stack 
        if (s.empty() || hist[s.top()] <= hist[i]) 
            s.push(i++); 
  
        // If this bar is lower than top of stack,  
        // then calculate area of rectangle with stack  
        // top as the smallest (or minimum height) bar.  
        // 'i' is 'right index' for the top and element  
        // before top in stack is 'left index' 
        else
        { 
            tp = s.top();  // store the top index 
            s.pop();  // pop the top 
  
            // Calculate the area with hist[tp] stack  
            // as smallest bar 
            area_with_top = hist[tp] * (s.empty() ? i :  
                                   i - s.top() - 1); 
  
            // update max area, if needed 
            if (max_area < area_with_top) 
                max_area = area_with_top; 
        } 
    } 
  
    // Now pop the remaining bars from stack and calculate 
    // area with every popped bar as the smallest bar 
    while (s.empty() == false) 
    { 
        tp = s.top(); 
        s.pop(); 
        area_with_top = hist[tp] * (s.empty() ? i :  
                                i - s.top() - 1); 
  
        if (max_area < area_with_top) 
            max_area = area_with_top; 
    } 
  
    return max_area; 
} 
  
// Driver program to test above function 
int main() 
{ 
    int hist[] = {6, 2, 5, 4, 5, 1, 6}; 
    int n = sizeof(hist)/sizeof(hist[0]); 
    cout << "Maximum area is " << getMaxArea(hist, n); 
    return 0; 
} 






// a linear time solution for stock span problem
#include <iostream>
#include <stack>
using namespace std;
  
// A stack based efficient method to calculate stock span values
void calculateSpan(int price[], int n, int S[])
{
   // Create a stack and push index of first element to it
   stack<int> st;

    st.push(0);
  
   // Span value of first element is always 1
   S[0] = 1;
  
   // Calculate span values for rest of the elements
   for (int i = 1; i < n; i++)
   {
      // Pop elements from stack while stack is not empty and top of
      // stack is smaller than price[i]
      while (!st.empty() && price[st.top()] <= price[i])
         st.pop();
  
      // If stack becomes empty, then price[i] is greater than all elements
      // on left of it, i.e., price[0], price[1],..price[i-1].  Else price[i]
      // is greater than elements after top of stack
      S[i] = (st.empty())? (i + 1) : (i - st.top());
  
      // Push this element to stack
      st.push(i);
   }
}




// C++ Program to find best buying and selling days 
#include <bits/stdc++.h> 
using namespace std; 
  
// solution structure 
class Interval { 
public: 
    int buy; 
    int sell; 
}; 
  
// This function finds the buy sell 
// schedule for maximum profit 
void stockBuySell(int price[], int n) 
{ 
    // Prices must be given for at least two days 
    if (n == 1) 
        return; 
  
    int count = 0; // count of solution pairs 
  
    // solution vector 
    Interval sol[n / 2 + 1]; 
  
    // Traverse through given price array 
    int i = 0; 
    while (i < n - 1) { 
        // Find Local Minima. Note that the limit is (n-2) as we are 
        // comparing present element to the next element. 
        while ((i < n - 1) && (price[i + 1] <= price[i])) 
            i++; 
  
        // If we reached the end, break 
        // as no further solution possible 
        if (i == n - 1) 
            break; 
  
        // Store the index of minima 
        sol[count].buy = i++; 
  
        // Find Local Maxima. Note that the limit is (n-1) as we are 
        // comparing to previous element 
        while ((i < n) && (price[i] >= price[i - 1])) 
            i++; 
  
        // Store the index of maxima 
        sol[count].sell = i - 1; 
  
        // Increment count of buy/sell pairs 
        count++; 
    } 
  
    // print solution 
    if (count == 0) 
        cout << "There is no day when buying"
             << " the stock will make profitn"; 
    else { 
        for (int i = 0; i < count; i++) 
            cout << "Buy on day: " << sol[i].buy 
                 << "\t Sell on day: " << sol[i].sell << endl; 
    } 
  
    return; 
} 
  
// Driver code 
int main() 
{ 
    // stock prices on consecutive days 
    int price[] = { 100, 180, 260, 310, 40, 535, 695 }; 
    int n = sizeof(price) / sizeof(price[0]); 
  
    // fucntion call 
    stockBuySell(price, n); 
  
    return 0; 
} 
  



// C++ program to find minimum breaks needed
// to break a string in dictionary words.
#include <bits/stdc++.h>
using namespace std;
  
const int ALPHABET_SIZE = 26;
  
// trie node
struct TrieNode {
    struct TrieNode* children[ALPHABET_SIZE];
  
    // isEndOfWord is true if the node 
    // represents end of a word
    bool isEndOfWord;
};
  
// Returns new trie node (initialized to NULLs)
struct TrieNode* getNode(void)
{
    struct TrieNode* pNode = new TrieNode;
  
    pNode->isEndOfWord = false;
  
    for (int i = 0; i < ALPHABET_SIZE; i++)
        pNode->children[i] = NULL;
  
    return pNode;
}
  
// If not present, inserts the key into the trie
// If the key is the prefix of trie node, just
// marks leaf node
void insert(struct TrieNode* root, string key)
{
    struct TrieNode* pCrawl = root;
  
    for (int i = 0; i < key.length(); i++) {
        int index = key[i] - 'a';
        if (!pCrawl->children[index])

        	  pCrawl->children[index] = getNode();
  
        pCrawl = pCrawl->children[index];
    }
  
    // mark last node as leaf
    pCrawl->isEndOfWord = true;
}
  
// function break the string into minimum cut
// such the every substring after breaking 
// in the dictionary.
void minWordBreak(struct TrieNode* root, 
          string key, int start, int* min_Break, 
                                 int level = 0)
{
    struct TrieNode* pCrawl = root;
  
    // base case, update minimum Break
    if (start == key.length()) {        
        *min_Break = min(*min_Break, level - 1);
        return;
    }
  
    // traverse given key(pattern)
    int minBreak = 0;   
    for (int i = start; i < key.length(); i++) {
        int index = key[i] - 'a';
        if (!pCrawl->children[index])
            return;
  
        // if we find a condition were we can 
        // move to the next word in a trie
        // dictionary
        if (pCrawl->children[index]->isEndOfWord)
            minWordBreak(root, key, i + 1,
                           min_Break, level + 1);
  
        pCrawl = pCrawl->children[index];
    }
}
  
// Driver program to test above functions
int main()
{
    string dictionary[] = { "Cat", "Mat",
   "Ca", "Ma", "at", "C", "Dog", "og", "Do" };
    int n = sizeof(dictionary) / sizeof(dictionary[0]);

        struct TrieNode* root = getNode();
  
    // Construct trie
    for (int i = 0; i < n; i++)
        insert(root, dictionary[i]);
    int min_Break = INT_MAX;
  
    minWordBreak(root, "CatMatat", 0, &min_Break, 0);
    cout << min_Break << endl;
    return 0;
}






// Suffix tree Implementation

// A simple C++ implementation of substring search using trie of suffixes
#include<iostream>
#include<list>
#define MAX_CHAR 256

using namespace std;
  
// A Suffix Trie (A Trie of all suffixes) Node
class SuffixTrieNode
{
private:
    SuffixTrieNode *children[MAX_CHAR];
    list<int> *indexes;
public:
    SuffixTrieNode() // Constructor
    {
        // Create an empty linked list for indexes of
        // suffixes starting from this node
        indexes = new list<int>;
  
        // Initialize all child pointers as NULL
        for (int i = 0; i < MAX_CHAR; i++)
          children[i] = NULL;
    }
  
    // A recursive function to insert a suffix of the txt
    // in subtree rooted with this node
    void insertSuffix(string suffix, int index);
  
    // A function to search a pattern in subtree rooted
    // with this node.The function returns pointer to a linked
    // list containing all indexes where pattern is present.
    // The returned indexes are indexes of last characters
    // of matched text.
    list<int>* search(string pat);
};
  
// A Trie of all suffixes
class SuffixTrie
{
private:
    SuffixTrieNode root;
public:
    // Constructor (Builds a trie of suffies of the given text)
    SuffixTrie(string txt)
    {
        // Consider all suffixes of given string and insert
        // them into the Suffix Trie using recursive function
        // insertSuffix() in SuffixTrieNode class
        for (int i = 0; i < txt.length(); i++)
            root.insertSuffix(txt.substr(i), i);
    }


       // Function to searches a pattern in this suffix trie.
    void search(string pat);
};
  
// A recursive function to insert a suffix of the txt in
// subtree rooted with this node
void SuffixTrieNode::insertSuffix(string s, int index)
{
    // Store index in linked list
    indexes->push_back(index);
  
    // If string has more characters
    if (s.length() > 0)
    {
        // Find the first character
        char cIndex = s.at(0);
  
        // If there is no edge for this character, add a new edge
        if (children[cIndex] == NULL)
            children[cIndex] = new SuffixTrieNode();
  
        // Recur for next suffix
        children[cIndex]->insertSuffix(s.substr(1), index+1);
    }
}
  
// A recursive function to search a pattern in subtree rooted with
// this node
list<int>* SuffixTrieNode::search(string s)
{
    // If all characters of pattern have been processed,
    if (s.length() == 0)
        return indexes;
  
    // if there is an edge from the current node of suffix trie,
    // follow the edge.
    if (children[s.at(0)] != NULL)
        return (children[s.at(0)])->search(s.substr(1));
  
    // If there is no edge, pattern doesn’t exist in text
    else return NULL;
}
  
/* Prints all occurrences of pat in the Suffix Trie S (built for text)*/
void SuffixTrie::search(string pat)
{
    // Let us call recursive search function for root of Trie.
    // We get a list of all indexes (where pat is present in text) in

     // variable 'result'
    list<int> *result = root.search(pat);
  
    // Check if the list of indexes is empty or not
    if (result == NULL)
        cout << "Pattern not found" << endl;
    else
    {
       list<int>::iterator i;
       int patLen = pat.length();
       for (i = result->begin(); i != result->end(); ++i)
         cout << "Pattern found at position " << *i - patLen<< endl;
    }
}
  
// driver program to test above functions
int main()
{
    // Let us build a suffix trie for text "geeksforgeeks.org"
    string txt = "geeksforgeeks.org";
    SuffixTrie S(txt);
  
    cout << "Search for 'ee'" << endl;
    S.search("ee");
  
    cout << "\nSearch for 'geek'" << endl;
    S.search("geek");
  
    cout << "\nSearch for 'quiz'" << endl;
    S.search("quiz");
  
    cout << "\nSearch for 'forgeeks'" << endl;
    S.search("forgeeks");
  
    return 0;
}







// C++ program to demonstrate auto-complete feature 
// using Trie data structure. 
#include<bits/stdc++.h> 
using namespace std; 
  
// Alphabet size (# of symbols) 
#define ALPHABET_SIZE (26) 
  
// Converts key current character into index 
// use only 'a' through 'z' and lower case 
#define CHAR_TO_INDEX(c) ((int)c - (int)'a') 
  
// trie node 
struct TrieNode 
{ 
    struct TrieNode *children[ALPHABET_SIZE]; 
  
    // isWordEnd is true if the node represents 
    // end of a word 
    bool isWordEnd; 
}; 
  
// Returns new trie node (initialized to NULLs) 
struct TrieNode *getNode(void) 
{ 
    struct TrieNode *pNode = new TrieNode; 
    pNode->isWordEnd = false; 
  
    for (int i = 0; i < ALPHABET_SIZE; i++) 
        pNode->children[i] = NULL; 
  
    return pNode; 
} 
  
// If not present, inserts key into trie.  If the 
// key is prefix of trie node, just marks leaf node 
void insert(struct TrieNode *root,  const string key) 
{ 
    struct TrieNode *pCrawl = root; 
  
    for (int level = 0; level < key.length(); level++) 
    { 
        int index = CHAR_TO_INDEX(key[level]); 
        if (!pCrawl->children[index]) 
            pCrawl->children[index] = getNode(); 
  
        pCrawl = pCrawl->children[index]; 
    } 
  
    // mark last node as leaf 
    pCrawl->isWordEnd = true; 
} 
  
// Returns true if key presents in trie, else false 
bool search(struct TrieNode *root, const string key) 
{ 
    int length = key.length(); 
    struct TrieNode *pCrawl = root; 
    for (int level = 0; level < length; level++) 
    { 
        int index = CHAR_TO_INDEX(key[level]); 
  
        if (!pCrawl->children[index]) 
            return false; 
  
        pCrawl = pCrawl->children[index]; 
    } 
  
    return (pCrawl != NULL && pCrawl->isWordEnd); 
} 
  
// Returns 0 if current node has a child 
// If all children are NULL, return 1. 
bool isLastNode(struct TrieNode* root) 
{ 
    for (int i = 0; i < ALPHABET_SIZE; i++) 
        if (root->children[i]) 
            return 0; 
    return 1; 
} 
  
// Recursive function to print auto-suggestions for given 
// node. 
void suggestionsRec(struct TrieNode* root, string currPrefix) 
{ 
    // found a string in Trie with the given prefix 
    if (root->isWordEnd) 
    { 
        cout << currPrefix; 
        cout << endl; 
    } 
  
    // All children struct node pointers are NULL 
    if (isLastNode(root)) 
        return; 
  
    for (int i = 0; i < ALPHABET_SIZE; i++) 
    { 
        if (root->children[i]) 
        { 
            // append current character to currPrefix string 
            currPrefix.push_back(97 + i); 
  
            // recur over the rest 
            suggestionsRec(root->children[i], currPrefix); 
        } 
    } 
} 
  
// print suggestions for given query prefix. 
int printAutoSuggestions(TrieNode* root, const string query) 
{ 
    struct TrieNode* pCrawl = root; 
  
    // Check if prefix is present and find the 
    // the node (of last level) with last character 
    // of given string. 
    int level; 
    int n = query.length(); 
    for (level = 0; level < n; level++) 
    { 
        int index = CHAR_TO_INDEX(query[level]); 
  
        // no string in the Trie has this prefix 
        if (!pCrawl->children[index]) 
            return 0; 
  
        pCrawl = pCrawl->children[index]; 
    } 
  
    // If prefix is present as a word. 
    bool isWord = (pCrawl->isWordEnd == true); 
  
    // If prefix is last node of tree (has no 
    // children) 
    bool isLast = isLastNode(pCrawl); 
  
    // If prefix is present as a word, but 
    // there is no subtree below the last 
    // matching node. 
    if (isWord && isLast) 
    { 
        cout << query << endl; 
        return -1; 
    } 
  
    // If there are are nodes below last 
    // matching character. 
    if (!isLast) 
    { 
        string prefix = query; 
        suggestionsRec(pCrawl, prefix); 
        return 1; 
    } 
} 
  
// Driver Code 
int main() 
{ 
    struct TrieNode* root = getNode(); 
    insert(root, "hello"); 
    insert(root, "dog"); 
    insert(root, "hell"); 
    insert(root, "cat"); 
    insert(root, "a"); 
    insert(root, "hel"); 
    insert(root, "help"); 
    insert(root, "helps"); 
    insert(root, "helping"); 
    int comp = printAutoSuggestions(root, "hel"); 
  
    if (comp == -1) 
        cout << "No other strings found with this prefix\n"; 
  
    else if (comp == 0) 
        cout << "No string found with this prefix\n"; 
  
    return 0; 
} 





// C++ program to find out maximum profit by 
// buying and selling a share atmost k times 
// given stock price of n days 
#include <climits> 
#include <iostream> 
using namespace std; 
  
// Function to find out maximum profit by buying 
// & selling a share atmost k times given stock 
// price of n days 
int maxProfit(int price[], int n, int k) 
{ 
    // table to store results of subproblems 
    // profit[t][i] stores maximum profit using 
    // atmost t transactions up to day i (including 
    // day i) 
    int profit[k + 1][n + 1]; 
  
    // For day 0, you can't earn money 
    // irrespective of how many times you trade 
    for (int i = 0; i <= k; i++) 
        profit[i][0] = 0; 
  
    // profit is 0 if we don't do any transation 
    // (i.e. k =0) 
    for (int j = 0; j <= n; j++) 
        profit[0][j] = 0; 
  
    // fill the table in bottom-up fashion 
    for (int i = 1; i <= k; i++) { 
        for (int j = 1; j < n; j++) { 
            int max_so_far = INT_MIN; 
  
            for (int m = 0; m < j; m++) 
                max_so_far = max(max_so_far, 
                                 price[j] - price[m] + profit[i - 1][m]); 
  
            profit[i][j] = max(profit[i][j - 1], max_so_far); 
        } 
    } 
  
    return profit[k][n - 1]; 
} 
  
// Driver code 
int main() 
{ 
    int k = 2; 
    int price[] = { 10, 22, 5, 75, 65, 80 }; 
    int n = sizeof(price) / sizeof(price[0]); 
  
    cout << "Maximum profit is: "
         << maxProfit(price, n, k); 
  
    return 0; 
} 



// Function to skip M nodes and then delete N nodes of the linked list.
void skipMdeleteN(struct Node  *head, int M, int N)
{
    struct Node *curr = head, *t;
    int count;
  
    // The main loop that traverses through the whole list


    while (curr)
    {
        // Skip M nodes
        for (count = 1; count<M && curr!= NULL; count++)
            curr = curr->next;
  
        // If we reached end of list, then return
        if (curr == NULL)
            return;
  
        // Start from next node and delete N nodes
        t = curr->next;
        for (count = 1; count<=N && t!= NULL; count++)
        {
            struct node *temp = t;
            t = t->next;
            free(temp);
        }
        curr->next = t; // Link the previous list with remaining nodes
  
        // Set current pointer for next iteration
        curr = t;
    }
}



// C++ implementation of the above approach 
#include <bits/stdc++.h> 
using namespace std; 
  
// Function to return minimum number of 
// operations to convert string A to B 
int minOperations(string s, string f) 
{ 
    unordered_map<string, int> vis; 
  
    int n; 
  
    n = s.length(); 
    int pos = 0; 
    for (int i = 0; i < s.length(); i++) { 
        if (s[i] == '_') { 
  
            // store the position of '_' 
            pos = i; 
            break; 
        } 
    } 
  
    // to store the generated string at every 
    // move and the position of '_' within it 
    queue<pair<string, int> > q; 
  
    q.push({ s, pos }); 
  
    // vis will store the minimum operations 
    // to reach that particular string 
    vis[s] = 0; 
  
    while (!q.empty()) { 
        string ss = q.front().first; 
        int pp = q.front().second; 
  
        // minimum moves to reach the string ss 
        int dist = vis[ss]; 
        q.pop(); 
  
        // try all 4 possible operations 
  
        // if '_' can be swapped with 
        // the character on it's left 
        if (pp > 0) { 
  
            // swap with the left character 
            swap(ss[pp], ss[pp - 1]); 
  
            // if the string is generated 
            // for the first time 
            if (!vis.count(ss)) { 
  
                // if generated string is 
                // the required string 
                if (ss == f) { 
                    return dist + 1; 
                    break; 
                } 
  
                // update the distance for the 
                // currently generated string 
                vis[ss] = dist + 1; 
                q.push({ ss, pp - 1 }); 
            } 
  
            // restore the string before it was 
            // swapped to check other cases 
            swap(ss[pp], ss[pp - 1]); 
        } 
  
        // swap '_' with the character 
        // on it's right this time 
        if (pp < n - 1) { 
            swap(ss[pp], ss[pp + 1]); 
            if (!vis.count(ss)) { 
                if (ss == f) { 
                    return dist + 1; 
                    break; 
                } 
                vis[ss] = dist + 1; 
                q.push({ ss, pp + 1 }); 
            } 
            swap(ss[pp], ss[pp + 1]); 
        } 
  
        // if '_' can be swapped 
        // with the character 'i+2' 
        if (pp > 1 && ss[pp - 1] != ss[pp - 2]) { 
            swap(ss[pp], ss[pp - 2]); 
            if (!vis.count(ss)) { 
                if (ss == f) { 
                    return dist + 1; 
                    break; 
                } 
                vis[ss] = dist + 1; 
                q.push({ ss, pp - 2 }); 
            } 
            swap(ss[pp], ss[pp - 2]); 
        } 
  
        // if '_' can be swapped 
        // with the character at 'i+2' 
        if (pp < n - 2 && ss[pp + 1] != ss[pp + 2]) { 
            swap(ss[pp], ss[pp + 2]); 
            if (!vis.count(ss)) { 
                if (ss == f) { 
                    return dist + 1; 
                    break; 
                } 
                vis[ss] = dist + 1; 
                q.push({ ss, pp + 2 }); 
            } 
            swap(ss[pp], ss[pp + 2]); 
        } 
    } 
} 
  
// Driver code 
int main() 
{ 
  
    string A = "aba_a"; 
    string B = "_baaa"; 
  
    cout << minOperations(A, B); 
  
    return 0; 
} 




#include <stdio.h> 
#include <stdlib.h> 
  
#define ARRAY_SIZE(arr) sizeof(arr)/sizeof(arr[0]) 
  
typedef struct node_t node_t; 
  
/* Binary tree node */
struct node_t 
{ 
    int data; 
    int lCount; 
  
    node_t* left; 
    node_t* right; 
}; 
  
/* Iterative insertion 
   Recursion is least preferred unless we gain something 
*/
node_t *insert_node(node_t *root, node_t* node) 
{ 
    /* A crawling pointer */
    node_t *pTraverse = root; 
    node_t *currentParent = root; 
  
    // Traverse till appropriate node 
    while(pTraverse) 
    { 
        currentParent = pTraverse; 
  
        if( node->data < pTraverse->data ) 
        { 
            /* We are branching to left subtree 
               increment node count */
            pTraverse->lCount++; 
            /* left subtree */
            pTraverse = pTraverse->left; 
        } 
        else
        { 
            /* right subtree */
            pTraverse = pTraverse->right; 
        } 
    } 
  
    /* If the tree is empty, make it as root node */
    if( !root ) 
    { 
        root = node; 
    } 
    else if( node->data < currentParent->data ) 
    { 
        /* Insert on left side */
        currentParent->left = node; 
    } 
    else
    { 
        /* Insert on right side */
        currentParent->right = node; 
    } 
  
    return root; 
} 
  
/* Elements are in an array. The function builds binary tree */
node_t* binary_search_tree(node_t *root, int keys[], int const size) 
{ 
    int iterator; 
    node_t *new_node = NULL; 
  
    for(iterator = 0; iterator < size; iterator++) 
    { 
        new_node = (node_t *)malloc( sizeof(node_t) ); 
  
        /* initialize */
        new_node->data   = keys[iterator]; 
        new_node->lCount = 0; 
        new_node->left   = NULL; 
        new_node->right  = NULL; 
  
        /* insert into BST */
        root = insert_node(root, new_node); 
    } 
  
    return root; 
} 
  
int k_smallest_element(node_t *root, int k) 
{ 
    int ret = -1; 
  
    if( root ) 
    { 
        /* A crawling pointer */
        node_t *pTraverse = root; 
  
        /* Go to k-th smallest */
        while(pTraverse) 
        { 
            if( (pTraverse->lCount + 1) == k ) 
            { 
                ret = pTraverse->data; 
                break; 
            } 
            else if( k > pTraverse->lCount ) 
            { 
                /*  There are less nodes on left subtree 
                    Go to right subtree */
                k = k - (pTraverse->lCount + 1); 
                pTraverse = pTraverse->right; 
            } 
            else
            { 
                /* The node is on left subtree */
                pTraverse = pTraverse->left; 
            } 
        } 
    } 
  
    return ret; 
} 
  
int main(void) 
{ 
    /* just add elements to test */
    /* NOTE: A sorted array results in skewed tree */
    int ele[] = { 20, 8, 22, 4, 12, 10, 14 }; 
    int i; 
    node_t* root = NULL; 
  
    /* Creating the tree given in the above diagram */
    root = binary_search_tree(root, ele, ARRAY_SIZE(ele)); 
  
    /*  It should print the sorted array */
    for(i = 1; i <= ARRAY_SIZE(ele); i++) 
    { 
        printf("\n kth smallest element for k = %d is %d", 
                 i, k_smallest_element(root, i)); 
    } 
  
    getchar(); 
    return 0; 
} 






// Program to show segment tree to demonstrate lazy 
// propagation 
#include <stdio.h> 
#include <math.h> 
#define MAX 1000 
  
// Ideally, we should not use global variables and large 
// constant-sized arrays, we have done it here for simplicity. 
int tree[MAX] = {0};  // To store segment tree 
int lazy[MAX] = {0};  // To store pending updates 
  
/*  si -> index of current node in segment tree 
    ss and se -> Starting and ending indexes of elements for 
                 which current nodes stores sum. 
    us and ue -> starting and ending indexes of update query 
    diff -> which we need to add in the range us to ue */
void updateRangeUtil(int si, int ss, int se, int us, 
                     int ue, int diff) 
{ 
    // If lazy value is non-zero for current node of segment 
    // tree, then there are some pending updates. So we need 
    // to make sure that the pending updates are done before 
    // making new updates. Because this value may be used by 
    // parent after recursive calls (See last line of this 
    // function) 
    if (lazy[si] != 0) 
    { 
        // Make pending updates using value stored in lazy 
        // nodes 
        tree[si] += (se-ss+1)*lazy[si]; 
  
        // checking if it is not leaf node because if 
        // it is leaf node then we cannot go further 
        if (ss != se) 
        { 
            // We can postpone updating children we don't 
            // need their new values now. 
            // Since we are not yet updating children of si, 
            // we need to set lazy flags for the children 
            lazy[si*2 + 1]   += lazy[si]; 
            lazy[si*2 + 2]   += lazy[si]; 
        } 
  
        // Set the lazy value for current node as 0 as it 
        // has been updated 
        lazy[si] = 0; 
    } 
  
    // out of range 
    if (ss>se || ss>ue || se<us) 
        return ; 
  
    // Current segment is fully in range 
    if (ss>=us && se<=ue) 
    { 
        // Add the difference to current node 
        tree[si] += (se-ss+1)*diff; 
  
        // same logic for checking leaf node or not 
        if (ss != se) 
        { 
            // This is where we store values in lazy nodes, 
            // rather than updating the segment tree itelf 
            // Since we don't need these updated values now 
            // we postpone updates by storing values in lazy[] 
            lazy[si*2 + 1]   += diff; 
            lazy[si*2 + 2]   += diff; 
        } 
        return; 
    } 
  
    // If not completely in rang, but overlaps, recur for 
    // children, 
    int mid = (ss+se)/2; 
    updateRangeUtil(si*2+1, ss, mid, us, ue, diff); 
    updateRangeUtil(si*2+2, mid+1, se, us, ue, diff); 
  
    // And use the result of children calls to update this 
    // node 
    tree[si] = tree[si*2+1] + tree[si*2+2]; 
} 
  
// Function to update a range of values in segment 
// tree 
/*  us and eu -> starting and ending indexes of update query 
    ue  -> ending index of update query 
    diff -> which we need to add in the range us to ue */
void updateRange(int n, int us, int ue, int diff) 
{ 
   updateRangeUtil(0, 0, n-1, us, ue, diff); 
} 
  
  
/*  A recursive function to get the sum of values in given 
    range of the array. The following are parameters for 
    this function. 
    si --> Index of current node in the segment tree. 
           Initially 0 is passed as root is always at' 
           index 0 
    ss & se  --> Starting and ending indexes of the 
                 segment represented by current node, 
                 i.e., tree[si] 
    qs & qe  --> Starting and ending indexes of query 
                 range */
int getSumUtil(int ss, int se, int qs, int qe, int si) 
{ 
    // If lazy flag is set for current node of segment tree, 
    // then there are some pending updates. So we need to 
    // make sure that the pending updates are done before 
    // processing the sub sum query 
    if (lazy[si] != 0) 
    { 
        // Make pending updates to this node. Note that this 
        // node represents sum of elements in arr[ss..se] and 
        // all these elements must be increased by lazy[si] 
        tree[si] += (se-ss+1)*lazy[si]; 
  
        // checking if it is not leaf node because if 
        // it is leaf node then we cannot go further 
        if (ss != se) 
        { 
            // Since we are not yet updating children os si, 
            // we need to set lazy values for the children 
            lazy[si*2+1] += lazy[si]; 
            lazy[si*2+2] += lazy[si]; 
        } 
  
        // unset the lazy value for current node as it has 
        // been updated 
        lazy[si] = 0; 
    } 
  
    // Out of range 
    if (ss>se || ss>qe || se<qs) 
        return 0; 
  
    // At this point we are sure that pending lazy updates 
    // are done for current node. So we can return value  
    // (same as it was for query in our previous post) 
  
    // If this segment lies in range 
    if (ss>=qs && se<=qe) 
        return tree[si]; 
  
    // If a part of this segment overlaps with the given 
    // range 
    int mid = (ss + se)/2; 
    return getSumUtil(ss, mid, qs, qe, 2*si+1) + 
           getSumUtil(mid+1, se, qs, qe, 2*si+2); 
} 
  
// Return sum of elements in range from index qs (quey 
// start) to qe (query end).  It mainly uses getSumUtil() 
int getSum(int n, int qs, int qe) 
{ 
    // Check for erroneous input values 
    if (qs < 0 || qe > n-1 || qs > qe) 
    { 
        printf("Invalid Input"); 
        return -1; 
    } 
  
    return getSumUtil(0, n-1, qs, qe, 0); 
} 
  
// A recursive function that constructs Segment Tree for 
//  array[ss..se]. si is index of current node in segment 
// tree st. 
void constructSTUtil(int arr[], int ss, int se, int si) 
{ 
    // out of range as ss can never be greater than se 
    if (ss > se) 
        return ; 
  
    // If there is one element in array, store it in 
    // current node of segment tree and return 
    if (ss == se) 
    { 
        tree[si] = arr[ss]; 
        return; 
    } 
  
    // If there are more than one elements, then recur 
    // for left and right subtrees and store the sum 
    // of values in this node 
    int mid = (ss + se)/2; 
    constructSTUtil(arr, ss, mid, si*2+1); 
    constructSTUtil(arr, mid+1, se, si*2+2); 
  
    tree[si] = tree[si*2 + 1] + tree[si*2 + 2]; 
} 
  
/* Function to construct segment tree from given array. 
   This function allocates memory for segment tree and 
   calls constructSTUtil() to fill the allocated memory */
void constructST(int arr[], int n) 
{ 
    // Fill the allocated memory st 
    constructSTUtil(arr, 0, n-1, 0); 
} 
  
  
// Driver program to test above functions 
int main() 
{ 
    int arr[] = {1, 3, 5, 7, 9, 11}; 
    int n = sizeof(arr)/sizeof(arr[0]); 
  
    // Build segment tree from given array 
    constructST(arr, n); 
  
    // Print sum of values in array from index 1 to 3 
    printf("Sum of values in given range = %d\n", 
           getSum(n, 1, 3)); 
  
    // Add 10 to all nodes at indexes from 1 to 5. 
    updateRange(n, 1, 5, 10); 
  
    // Find sum after the value is updated 
    printf("Updated sum of values in given range = %d\n", 
            getSum( n, 1, 3)); 
  
    return 0; 
} 





// C++ program to find maximum revenue by placing 
// billboard on the highway with given constarints. 
#include<bits/stdc++.h> 
using namespace std; 
  
int maxRevenue(int m, int x[], int revenue[], int n, 
                                              int t) 
{ 
    // Array to store maximum revenue at each miles. 
    int maxRev[m+1]; 
    memset(maxRev, 0, sizeof(maxRev)); 
  
    // actual minimum distance between 2 billboards. 
    int nxtbb = 0; 
    for (int i = 1; i <= m; i++) 
    { 
        // check if all billboards are already placed. 
        if (nxtbb < n) 
        { 
            // check if we have billboard for that particular 
            // mile. If not, copy the previous maximum revenue. 
            if (x[nxtbb] != i) 
                maxRev[i] = maxRev[i-1]; 
  
            // we do have billboard for this mile. 
            else
            { 
                // We have 2 options, we either take current  
                // or we ignore current billboard. 
  
                // If current position is less than or equal to 
                // t, then we can have only one billboard. 
                if (i <= t) 
                    maxRev[i] = max(maxRev[i-1], revenue[nxtbb]); 
  
                // Else we may have to remove previously placed 
                // billboard 
                else
                    maxRev[i] = max(maxRev[i-t-1]+revenue[nxtbb], 
                                                  maxRev[i-1]); 
  
                nxtbb++; 
            } 
        } 
        else
            maxRev[i] = maxRev[i - 1]; 
    } 
  
    return maxRev[m]; 
} 
  
// Driven Program 
int main() 
{ 
    int m = 20; 
    int x[] = {6, 7, 12, 13, 14}; 
    int revenue[] = {5, 6, 5, 3, 1}; 
    int n = sizeof(x)/sizeof(x[0]); 
    int t = 5; 
    cout << maxRevenue(m, x, revenue, n, t) << endl; 
    return 0; 
} 







// Delete duplicates from string and also do lexicographical order


// public String removeDuplicateLetters(String s) {
//     int[] counters = new int[26];
//     boolean[] visited = new boolean[26];
//     Stack<Character> stack = new Stack<>();

//     for(char c : s.toCharArray())
//       counters[c - 'a']++;

//     // use combination of visited and monotonous increasing stack
//     for(char c : s.toCharArray()) {
//       if(visited[c - 'a']) {
//         counters[c - 'a']--;
//         continue;
//       }

//       while(!stack.isEmpty() &&
//             stack.peek() > c &&
//             counters[stack.peek() - 'a'] > 1) {
//         counters[stack.peek() - 'a']--;
//         visited[stack.peek() - 'a'] = false;
//         stack.pop();
//       }

//       stack.push(c);
//       visited[c - 'a'] = true;
//     }

//     char[] result = new char[stack.size()];
//     for(int i = stack.size() - 1; i >= 0; i--)
//       result[i] = stack.pop();

//     return new String(result);
//   }



#include <iostream>
using namespace std;
int det[5][5];
int mat[13][5];

void detonate(int r)
{
    for(int i=r;i>r-5 && i>=0;i--)
    {
        for(int j=0;j<5;j++)
        {
            det[r-i][j]=0;
            if(mat[i][j]==2)
            {
                det[r-i][j]=2;
                mat[i][j]=0;
            }
        }
    }
}

void undet(int r)
{
    for(int i=r;i>r-5 && i>=0;i--)
        for(int j=0;j<5;j++)
        {
            if( det[r-i][j]==2)
                mat[i][j]=2;
        }
}
void func(int n,int pos,int c,int &max)
{
    if(pos>4||pos<0||c<0)
        return;

    if(mat[n][pos]==2)
        c-=1;
    else if(mat[n][pos]==1)
        c+=1;

    if(n==0)
    {
        if(c>max)
            max=c;
        return;
    }
    else
    {
        func(n-1,pos+1,c,max);
        func(n-1,pos-1,c,max);
        func(n-1,pos,c,max);
    }
}
int main()
{
    int t;
    cin>>t;
    int count=1;
    while(t--)
    {
        int n;
        cin>>n;
        for(int i=0;i<n;i++)
            for(int j=0;j<5;j++)
                cin>>mat[i][j];
        int max=-1,c;
        for(int j=0;j<5;j++)
            mat[n][j]=0;
        mat[n][2]=3;
        for(int j=n;j>=5;j--)
        {
            c=-1;
            detonate(j-1);
            func(n,2,0,c);
            if(c>max)
                max=c;
            undet(j-1);
        }
        if(max<0)
            max=-1;
        cout<<"#"<<count<<" "<<max<<endl;
        count++;
    }
}



// Wormhole.cpp : Defines the entry point for the console application.
//
 
#include <iostream>
using namespace std;
class point{
public:
    int x,y;
};
int find_min(int a[][100],int dist[],int vis[],int n)
{
    int min1=999999,index;
    for(int i=0;i<=n;i++)
    {
        if(vis[i]==0 && dist[i]<=min1)
            min1=dist[i],index=i;
    }
    return index;
}
int dijkstra(int a[][100],int n,int vis[])
{
    //vis[src]=1;
    int dist[100];
    //vis[0]=1;
    for(int i=0;i<=n;i++)
        dist[i]=999999;
    dist[0]=0;
    for(int i=0;i<=n;i++)
    {
        int temp=find_min(a,dist,vis,n);
        vis[temp]=1;
        for(int j=0;j<=n;j++)
        {
            if(!vis[j] && a[temp][j] && dist[temp]!=999999 && dist[j]>a[temp][j]+dist[temp])
                dist[j]=a[temp][j]+dist[temp];
        }
    //  cout<<dist[n]<<" ";
    }
    return dist[n];
}
int main()
{
    int t;
    cin>>t;
    for(int i=1;i<=t;i++)
    {
        int n;
        cin>>n;
        point p[100];
        int w[100];
        cin>>p[0].x>>p[0].y;
        cin>>p[2*n+1].x>>p[2*n+1].y;
        int j=1;
        for(int i=1;i<=2*n;i+=2)
        {
            cin>>p[i].x>>p[i].y>>p[i+1].x>>p[i+1].y>>w[j];
            j++;
        }
        int a[100][100]={{0}};
        for(int i=0;i<=2*n+1;i++)
        {
            for(int j=i;j<=2*n+1;j++)
            {
                a[i][j]=abs(p[i].x-p[j].x)+abs(p[i].y-p[j].y);
            }
        }
        int k=1;
        for(int i=1;i<=2*n;i+=2)
        {
            a[i][i+1]=w[k];
            a[i+1][i]=w[k];
            k++;
        }
        if(n==0)
        {
            cout<<abs(p[1].x-p[0].x)+abs(p[1].y-p[0].y)<<endl;
        }
        else
        {
            int vis[100]={0};
            cout<<dijkstra(a,2*n+1,vis)<<endl;
        }
    }
    return 0;
}
 














// C++ program to simplify algebraic string 
#include <bits/stdc++.h> 
using namespace std; 
  
// Function to simplify the string 
char* simplify(string str) 
{ 
    int len = str.length(); 
  
    // resultant string of max length equal 
    // to length of input string 
    char* res = new char(len); 
    int index = 0, i = 0; 
  
    // create empty stack 
    stack<int> s; 
    s.push(0); 
  
    while (i < len) { 
        if (str[i] == '+') { 
  
            // If top is 1, flip the operator 
            if (s.top() == 1) 
                res[index++] = '-'; 
  
            // If top is 0, append the same operator 
            if (s.top() == 0) 
                res[index++] = '+'; 
  
        } else if (str[i] == '-') { 
            if (s.top() == 1) 
                res[index++] = '+'; 
            else if (s.top() == 0) 
                res[index++] = '-'; 
        } else if (str[i] == '(' && i > 0) { 
            if (str[i - 1] == '-') { 
  
                // x is opposite to the top of stack 
                int x = (s.top() == 1) ? 0 : 1; 
                s.push(x); 
            } 
  
            // push value equal to top of the stack 
            else if (str[i - 1] == '+') 
                s.push(s.top()); 
        } 
  
        // If closing parentheses pop the stack once 
        else if (str[i] == ')') 
            s.pop(); 
  
        // copy the character to the result 
        else
            res[index++] = str[i]; 
        i++; 
    } 
    return res; 
} 
  
// Driver program 
int main() 
{ 
    string s1 = "a-(b+c)"; 
    string s2 = "a-(b-c-(d+e))-f"; 
    cout << simplify(s1) << endl; 
    cout << simplify(s2) << endl; 
    return 0; 
}



// C++ program to find sum of all minimum and maximum 
// elements Of Sub-array Size k. 
#include<bits/stdc++.h> 
using namespace std; 
  
// Returns sum of min and max element of all subarrays 
// of size k 
int SumOfKsubArray(int arr[] , int n , int k) 
{ 
    int sum = 0;  // Initialize result 
  
    // The queue will store indexes of useful elements 
    // in every window 
    // In deque 'G' we maintain decreasing order of 
    // values from front to rear 
    // In deque 'S' we  maintain increasing order of 
    // values from front to rear 
    deque< int > S(k), G(k); 
  
    // Process first window of size K 
    int i = 0; 
    for (i = 0; i < k; i++) 
    { 
        // Remove all previous greater elements 
        // that are useless. 
        while ( (!S.empty()) && arr[S.back()] >= arr[i]) 
            S.pop_back(); // Remove from rear 
  
        // Remove all previous smaller that are elements 
        // are useless. 
        while ( (!G.empty()) && arr[G.back()] <= arr[i]) 
            G.pop_back(); // Remove from rear 
  
        // Add current element at rear of both deque 
        G.push_back(i); 
        S.push_back(i); 
    } 
  
    // Process rest of the Array elements 
    for (  ; i < n; i++ ) 
    { 
        // Element at the front of the deque 'G' & 'S' 
        // is the largest and smallest 
        // element of previous window respectively 
        sum += arr[S.front()] + arr[G.front()]; 
  
        // Remove all elements which are out of this 
        // window 
        while ( !S.empty() && S.front() <= i - k) 
            S.pop_front(); 
        while ( !G.empty() && G.front() <= i - k) 
            G.pop_front(); 
  
        // remove all previous greater element that are 
        // useless 
        while ( (!S.empty()) && arr[S.back()] >= arr[i]) 
            S.pop_back(); // Remove from rear 
  
        // remove all previous smaller that are elements 
        // are useless 
        while ( (!G.empty()) && arr[G.back()] <= arr[i]) 
            G.pop_back(); // Remove from rear 
  
        // Add current element at rear of both deque 
        G.push_back(i); 
        S.push_back(i); 
    } 
  
    // Sum of minimum and maximum element of last window 
    sum += arr[S.front()] + arr[G.front()]; 
  
    return sum; 
} 
  
// Driver program to test above functions 
int main() 
{ 
    int arr[] = {2, 5, -1, 7, -3, -1, -2} ; 
    int n = sizeof(arr)/sizeof(arr[0]); 
    int k = 3; 
    cout << SumOfKsubArray(arr, n, k) ; 
    return 0; 
} 