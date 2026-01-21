dataset_name="G2";
data=readtable(strcat("/home/zkl/Desktop/9.3/",dataset_name,"/step1/act_gene_count.csv"),'VariableNamingRule','preserve');
edges=readtable(strcat("/home/zkl/Desktop/9.3/",dataset_name,"/step1/act_edge_name.xlsx"),"NumHeaderLines",1);

boxsize = 0.1;

[n_spot, ~] = size(data);
[n_edge, ~] = size(edges);
fprintf("%d\n",n_spot);
fprintf("%d\n",n_edge);

res = zeros(1,n_spot);
for i=1:n_edge
    fprintf("%d\n",i);
    gx = data.(edges(i,:).Var2{1,1});
    gy = data.(edges(i,:).Var3{1,1});
    res(i,:)=[csnedge(gx',gy',boxsize)];
end

csvwrite(strcat("/home/zkl/Desktop/9.3/",dataset_name,"/step2/act.csv"), res)
