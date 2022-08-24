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
%--------------------------------读取d_node_pos和d_node_reorder
d_dof = csvread('../../matrix/d_dof.csv');
d_node_pos = cell(4,1);
d_freenodesx = cell(4,1);
d_freenodesy = cell(4,1);
d_nodeid = cell(4,1);
d_node_reorder = cell(4,1);

for d = 1:4
    d_freenodesx{d} = zeros(d_dof(d),1);
    d_freenodesy{d} = zeros(d_dof(d),1);
    d_nodeid{d} = zeros(d_dof(d),1);
    d_node_reorder{d} = zeros(d_dof(d),1);
    str = sprintf('../../matrix/d_node_pos[%d].csv', d-1);
    d_node_pos{d} =  csvread(str) + 1;
    str = sprintf('../../matrix/d_nodeid[%d].csv', d-1);
    d_nodeid{d} =  csvread(str) + 1;
    str = sprintf('../../matrix/d_node_reorder[%d].csv', d-1);
    d_node_reorder{d} =  csvread(str) + 1;
    for n = 1:length(d_node_pos{d})
        if(d_node_pos{d}(n) <= d_dof(d))
            globalid = d_nodeid{d}(n);
            reorderid = d_node_pos{d}(n);
        end
    end
end
% %------------------------------------读取C++求解得到的S
% CS = cell(4,1);
% for i = 1:4
%     str = sprintf('../../matrix/S_domain[%d].csv',i-1);
%     CS{i} = csvread(str);
%     CS{i} = CS{i}(1:size(CS{i},1),1:size(CS{i},1));
% end
%------------------------------------绘制最终结果
At_real = csvread('../../matrix/At_real.csv');
figure(1);
plot(At_real,'ob','MarkerSize',3);
hold on;
%------------------------------------求解过程
deltaS = cell(4,1);
deltaF = cell(4,1);
res = cell(4,1);
errorSid = cell(4,1);
errorSreorder = cell(4,1);
errorFid = cell(4,1);
errorFreorder = cell(4,1);
At = zeros(length(Coor), 1);
At_old = zeros(length(Coor),1);
% %----------------------------------读取单元分区并绘制
% figure;
% DomainElement1 = cell(4,1);
% Domainx = cell(4,1);
% Domainy = cell(4,1);
% for i = 1:4
%     DomainElement1{i} = find(ePartTable(:, i));
%     Domainx{i} = zeros(length(DomainElement1{i}),3);
%     Domainy{i} = zeros(length(DomainElement1{i}),3);
%     for j = 1:3
%         Domainx{i}(:,j) = Coor(TriElement(DomainElement1{i},j),1);
%         Domainy{i}(:,j) = Coor(TriElement(DomainElement1{i},j),2);
%     end
% end
% patch(Domainx{1}',Domainy{1}','red','FaceAlpha',.3);
% hold on;
% patch(Domainx{2}',Domainy{2}','green','FaceAlpha',.3);
% hold on;
% patch(Domainx{3}',Domainy{3}','blue','FaceAlpha',.3);
% hold on;
% patch(Domainx{4}',Domainy{4}','yellow','FaceAlpha',.3);
% hold on;
% axis equal;
for iter = 1:200
    %------------------------------------读取C++求解得到的F
    CF = cell(4,1);
    for i = 1:4
        str = sprintf('../../matrix/F_iter[%d]_domain[%d].csv',iter-1,i-1);
        CF{i} = csvread(str);
    end
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
%                         if(d_n1 == 808 && d_n2 == 808)
%                             fprintf("d_n1: 808, globalnode: %d, x: %f, y: %f, eleid: %d, mu: %.12f, r: %.12f, Ce: %.12f, Se: %f\n", n1, Coor(n1,1), Coor(n1,2), eleID,mu(eleID),ydot(eleID),(R(eleID,row)*R(eleID,col)+Q(eleID,row)*Q(eleID,col))/4/AREA(eleID),CE(row,col));
%                         end
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
%         %--------------------------------------对比C++计算的S和MATLAB计算的S，对比错误的节点位置
%         deltaS{d} = S-CS{d};
%         deltaS{d}(deltaS{d} < 1e-5) = 0;
%         [errorSid{d}, errorSreorder{d}] = find(deltaS{d} ~= 0);
%         errorSid{d} = d_node_reorder{d}(errorSid{d});
% %         plot(Coor(d_nodeid{d}(errorrow{d}),1), Coor(d_nodeid{d}(errorrow{d}),2), '.r');
% %         hold on;
%         %--------------------------------------对比C++计算的S和MATLAB计算的F，对比错误的节点位置
%         deltaF{d} = F-CF{d};
%         deltaF{d}(deltaF{d} < 1e-5) = 0;
%         errorFreorder{d} = find(deltaF{d} ~= 0);
%         errorFid{d} = d_node_reorder{d}(errorFid{d});
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
    %---------------------------------对比求解结果
    str = sprintf('step: %d', iter);
    title(str);
    h1 = plot(At,'.r');
    str = sprintf('../../matrix/At_step%d.csv', iter-1);
    At_step = csvread(str);
    h2 = plot(At_step, 'og','MarkerSize',4);
     pause(0.5);
    delete(h1);
    delete(h2);
    %---------------------------------误差分析
    error = norm(At-At_old)/norm(At);
    fprintf('step: %d, error: %d\n', iter, error);
    if(error < 1e-5)
        return;
    else
        At_old = At;
    end
end