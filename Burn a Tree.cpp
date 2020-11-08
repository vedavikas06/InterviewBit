/**
 * Definition for binary tree
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
// Output format -> {present(1)/not-present(-1), {height of the sub tree, distance to the given leaf node having value B } }
 pair<int,pair<int,int>> burn_tree(TreeNode* A, int b, int &maxx){
    if(A==NULL){
        return {-1,{0,0}};
    }
    if(A->val==b){
        return {1,{1,1}};
    }
    pair<int,pair<int,int>> v1 = burn_tree(A->left, b, maxx);
    pair<int,pair<int,int>> v2 = burn_tree(A->right, b, maxx);
    
    if(v1.first==-1 && v2.first==-1){
        return {-1, {max(v1.second.first, v2.second.first)+1, 0}};
    }else if(v1.first!=-1 && v2.first==-1){
        maxx = max(maxx, v1.second.second+v2.second.first);
        return {1, {max(v1.second.first, v2.second.first)+1, v1.second.second+1}};
    }else{
        maxx = max(maxx, v2.second.second+v1.second.first);
        return {1, {max(v1.second.first, v2.second.first)+1, v2.second.second+1}};
    }
}

int Solution::solve(TreeNode* A, int B) {
    int maxx = 0;
    burn_tree(A,B,maxx);
    return maxx;
}

