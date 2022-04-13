% Schwarz区域分解
% By QzLancer
% 2022/4/13

%---------------------------------读取文件
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('model3656.mphtxt');
[fileID] = fopen('ContactLinear_7172.metis.epart.4');
eDomain = fscanf(fileID,'%d\n');
[fileID] = fopen('ContactLinear_7172.metis.npart.4');
nDomain = fscanf(fileID,'%d\n');
%----------------------------------计算单元的基本参数
num_elements = length(TriElement);
P = zeros(num_elements,3);
Q = zeros(num_elements,3);
R = zeros(num_elements,3);

X = Coor(:,1);
Y = Coor(:,2);
XL = X(TriElement);
YL = Y(TriElement);

Q(:,1) = YL(:,2) - YL(:,3);
Q(:,2) = YL(:,3) - YL(:,1);
Q(:,3) = YL(:,1) - YL(:,2);

R(:,1) = XL(:,3) - XL(:,2);
R(:,2) = XL(:,1) - XL(:,3);
R(:,3) = XL(:,2) - XL(:,1);

P(:,1) = XL(:,2).*YL(:,3) - XL(:,3).*YL(:,2);
P(:,2) = XL(:,3).*YL(:,1) - XL(:,1).*YL(:,3);
P(:,3) = XL(:,1).*YL(:,2) - XL(:,2).*YL(:,1);

AREA = 0.5 * (Q(:,1).*R(:,2) - Q(:,2).*R(:,1));%三角形面积

J = zeros(num_elements,1);%计算电流密度矩阵
coildomain = find(TriEntity == 5);%寻找线圈区域的单元
J(coildomain) = 8e6;%设置线圈区域的电流密度，其他其余为0

mu0 = 4*pi*1e-7;
mu = mu0*ones(num_elements,1);%保存每一个单元的磁导率，初始化为空气磁导率

ydot = zeros(num_elements,1);
for row=1:num_elements
    %两个点在坐标轴上,注意P59页公式是错误的，应当为x
    if XL(row,1)+XL(row,2)<1e-10 || XL(row,2)+XL(row,3)<1e-10 || XL(row,1)+XL(row,3)<1e-10
        ydot(row) = mean(XL(row,:));
    else
        ydot(row) = 1.5/(1/(XL(row,1)+XL(row,2))+1/(XL(row,1)+XL(row,3))+1/(XL(row,2)+XL(row,3)));
    end
end
%----------------------------------读取ePartTable和boundaryTable
ePartTable = csvread('../../matrix/eparttable.csv');
ePartTable = ePartTable(:,1:(length(ePartTable)-1))';
boundaryTable = csvread('../../matrix/boundarytable.csv');
boundaryTable = boundaryTable(:,1:length(boundaryTable))';
DomainElement = cell(4,1);
for d = 1:4
    DomainElement{d} = find(ePartTable(:, d) ~= -1);
end
DomainNode = cell(4,1); 
for d = 1:4
   DomainNode{d} = find(boundaryTable(:,d)); 
end
%--------------------------------读取nPartTable
nPartTable =  csvread('../../matrix/nparttable.csv');
nPartTable = nPartTable(:,1:length(Coor))';
%--------------------------------读取d_node_pos
d_dof = csvread('../../matrix/d_dof.csv');
d_node_pos = cell(4,1);
d_freenodesx = cell(4,1);
d_freenodesy = cell(4,1);
d_nodeid = cell(4,1);

for d = 1:4
    d_freenodesx{d} = zeros(d_dof(d),1);
    d_freenodesy{d} = zeros(d_dof(d),1);
    d_nodeid{d} = zeros(d_dof(d),1);
    str = sprintf('../../matrix/d_node_pos[%d].csv', d-1);
    d_node_pos{d} =  csvread(str) + 1;
    str = sprintf('../../matrix/d_nodeid[%d].csv', d-1);
    d_nodeid{d} =  csvread(str) + 1;
    for n = 1:length(d_node_pos{d})
        if(d_node_pos{d}(n) <= d_dof(d))
            globalid = d_nodeid{d}(n);
            reorderid = d_node_pos{d}(n);
        end
    end
end
%------------------------------------绘制真解
At_real = csvread('../../matrix/At_real.csv');
figure(5);
plot(At_real,'ob','MarkerSize',3);
hold on;
%------------------------------------求解过程
res = cell(4,1);
At = zeros(length(Coor), 1);
At_old = zeros(length(Coor),1);
for iter = 1:300
    %--------------------------------求解各个子域
    for d = 1:4
        res{d} = zeros(d_dof(d),1);
        CE = zeros(3,3);
        S = zeros(d_dof(d), d_dof(d));
        F = zeros(d_dof(d), 1);
        for e = 1:length(DomainElement{d})
            eleID = DomainElement{d}(e);
            for row = 1:3
                n1 = TriElement(eleID, row);
                d_n1 = nPartTable(n1,d) + 1;
                if d_node_pos{d}(d_n1) <= d_dof(d)
                    for col = 1:3
                        n2 = TriElement(eleID, col);
                        d_n2 = nPartTable(n2,d) + 1;
                        CE(row,col) = (R(eleID,row)*R(eleID,col)+Q(eleID,row)*Q(eleID,col))/4/AREA(eleID)/mu(eleID)/ydot(eleID);
                        if d_node_pos{d}(d_n2) <= d_dof(d)
                            S(d_node_pos{d}(d_n1), d_node_pos{d}(d_n2)) = S(d_node_pos{d}(d_n1), d_node_pos{d}(d_n2)) + CE(row, col);
                        else
                            F(d_node_pos{d}(d_n1)) = F(d_node_pos{d}(d_n1)) - CE(row, col)*At(n2);
                        end
                    end
                    F(d_node_pos{d}(d_n1)) = F(d_node_pos{d}(d_n1)) + J(eleID)*AREA(eleID)/3;
                end
            end
        end
        res{d} = S\F;
    end
    %--------------------------------子域结果整合到全局
    for d = 1:4
        for d_n = 1:length(d_nodeid{d})
            g_n = d_nodeid{d}(d_n);
            pos_n = d_node_pos{d}(d_n);
            if(pos_n <= d_dof(d))
%                 if At(g_n) ~= 0
%                     fprintf("d_n: %d, g_n: %d, pos_n: %d\n", d_n, g_n, pos_n);
%                 end
                At(g_n) = res{d}(pos_n);
            end
        end
    end
    str = sprintf('step: %d', iter);
    title(str);
    h = plot(At,'.r');
    pause(0.5);
    delete(h);
    error = norm(At-At_old)/norm(At);
    fprintf('step: %d, error: %d\n', iter, error);
    if(error < 1e-5)
        return;
    else
        At_old = At;
    end
end
h = plot(At,'.r');