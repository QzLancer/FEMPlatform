% metisֻ��ʵ�ַ��ص���������ֽ⣬ͨ��metis��������ݸ�ʽ���õ��ص����������ֽ�
% 4������
% metis�����ṩ������Ϣ��1����Ԫ���ڷ��ص������� 2���ڵ���������
% 1������ÿ���ڵ��ܱߵĵ�Ԫ���
% 2�������ڵ����ڵ����򣬽��ڵ���Χ�ĵ�Ԫ����Ǹýڵ����������ţ��õ���Ԫ�����ص�������
% 3�����������е����е�Ԫ�������Ԫ�Ľڵ���������͵�Ԫ��������һ����Ϊ�߽�ڵ�
% ����ȫ���ڵ㣬���������ýڵ��ȫ����Ԫ�����ɵõ��ص��汾������ֽ�
% By QzLancer
% 2022/3/28

%---------------------------------��ȡ�ļ�
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('model3656.mphtxt');
[fileID] = fopen('ContactLinear_7172.metis.epart.4');
eDomain = fscanf(fileID,'%d\n');
[fileID] = fopen('ContactLinear_7172.metis.npart.4');
nDomain = fscanf(fileID,'%d\n');
%--------------------------------��ȡ�ڵ���Χ��ȫ����Ԫ
EleSizeinNode = zeros(length(Coor),1);
GlobalEleIDinNode = zeros(length(Coor), 10);
for e = 1 : length(TriElement)
    for n = 1 :3
        nodeID = TriElement(e,n);
        EleSizeinNode(nodeID) = EleSizeinNode(nodeID) + 1;
        GlobalEleIDinNode(nodeID, EleSizeinNode(nodeID)) = e;
    end
end
%--------------------------------�����ڵ����ڵ����򣬽���ص�Ԫ����epartTable��
ePartTable = zeros(length(TriElement), 4);
for n = 1 : length(Coor)
    domain = nDomain(n) + 1;
    for e = 1 : EleSizeinNode(n)
        GlobalEleID = GlobalEleIDinNode(n, e);
        ePartTable(GlobalEleID, domain) = 1;
    end
end
%----------------------------------��ȡ��Ԫ����������
DomainElement = cell(4,1);
Domainx = cell(4,1);
Domainy = cell(4,1);
for i = 1:4
    DomainElement{i} = find(ePartTable(:, i));
    Domainx{i} = zeros(length(DomainElement{i}),3);
    Domainy{i} = zeros(length(DomainElement{i}),3);
    for j = 1:3
        Domainx{i}(:,j) = Coor(TriElement(DomainElement{i},j),1);
        Domainy{i}(:,j) = Coor(TriElement(DomainElement{i},j),2);
    end
end
patch(Domainx{1}',Domainy{1}','red','FaceAlpha',.3);
hold on;
patch(Domainx{2}',Domainy{2}','green','FaceAlpha',.3);
hold on;
patch(Domainx{3}',Domainy{3}','blue','FaceAlpha',.3);
hold on;
patch(Domainx{4}',Domainy{4}','yellow','FaceAlpha',.3);
hold on;
axis equal;
%----------------------------------����ÿ����Ԫ�ı߽�ڵ㣬������boundaryTable��
boundaryTable = zeros(length(Coor), 4);
for i = 1:4
    for e = 1:length(DomainElement{i})
        eleID = DomainElement{i}(e);
        for j = 1:3
            nodeID = TriElement(eleID,j);
            if nDomain(nodeID)+1 ~= i
                boundaryTable(nodeID, i) = 1;
            end
        end
    end
end
%----------------------------------����boundaryTable
DomainNode = cell(4,1); 
for i = 1:4
   DomainNode{i} = find(boundaryTable(:,i)); 
end
plot(Coor(DomainNode{1},1), Coor(DomainNode{1},2), '.r');
plot(Coor(DomainNode{2},1), Coor(DomainNode{2},2), '.g');
plot(Coor(DomainNode{3},1), Coor(DomainNode{3},2), '.b');
plot(Coor(DomainNode{4},1), Coor(DomainNode{4},2), '.y');
%----------------------------------��֤C++���ɵ�eparttable��boundaryTable�Ƿ���ȷ
ePartTable1 = csvread('../../matrix/eparttable.csv');
ePartTable1 = ePartTable1(:,1:length(ePartTable))';
boundaryTable1 = csvread('../../matrix/boundarytable.csv');
boundaryTable1 = boundaryTable1(:,1:length(boundaryTable))';
DomainElement = cell(4,1);
Domainx = cell(4,1);
Domainy = cell(4,1);
for i = 1:4
    DomainElement{i} = find(ePartTable1(:, i) ~= -1);
    Domainx{i} = zeros(length(DomainElement{i}),3);
    Domainy{i} = zeros(length(DomainElement{i}),3);
    for j = 1:3
        Domainx{i}(:,j) = Coor(TriElement(DomainElement{i},j),1);
        Domainy{i}(:,j) = Coor(TriElement(DomainElement{i},j),2);
    end
end
figure(2);
patch(Domainx{1}',Domainy{1}','red','FaceAlpha',.3);
hold on;
patch(Domainx{2}',Domainy{2}','green','FaceAlpha',.3);
hold on;
patch(Domainx{3}',Domainy{3}','blue','FaceAlpha',.3);
hold on;
patch(Domainx{4}',Domainy{4}','yellow','FaceAlpha',.3);
hold on;
axis equal;
DomainNode = cell(4,1); 
for i = 1:4
   DomainNode{i} = find(boundaryTable1(:,i)); 
end
plot(Coor(DomainNode{1},1), Coor(DomainNode{1},2), '.r');
plot(Coor(DomainNode{2},1), Coor(DomainNode{2},2), '.g');
plot(Coor(DomainNode{3},1), Coor(DomainNode{3},2), '.b');
plot(Coor(DomainNode{4},1), Coor(DomainNode{4},2), '.y');
%--------------------------------��֤C++���ɵ�npartTable�Ƿ���ȷ
nPartTable =  csvread('../../matrix/nparttable.csv');
nPartTable = nPartTable(:,1:length(Coor))';
DomainNode = cell(4,1);
Domainx = cell(4,1);
Domainy = cell(4,1);
for i = 1:4
    DomainElement{i} = find(ePartTable1(:, i) ~= -1);
    Domainx{i} = zeros(length(DomainElement{i}),3);
    Domainy{i} = zeros(length(DomainElement{i}),3);
    for j = 1:3
        Domainx{i}(:,j) = Coor(TriElement(DomainElement{i},j),1);
        Domainy{i}(:,j) = Coor(TriElement(DomainElement{i},j),2);
    end
end
figure(3);
patch(Domainx{1}',Domainy{1}','red','FaceAlpha',.3);
hold on;
patch(Domainx{2}',Domainy{2}','green','FaceAlpha',.3);
hold on;
patch(Domainx{3}',Domainy{3}','blue','FaceAlpha',.3);
hold on;
patch(Domainx{4}',Domainy{4}','yellow','FaceAlpha',.3);
hold on;
axis equal;
for i = 1:4
    DomainNode{i} = find(nPartTable(:, i) ~= -1);
end
plot(Coor(DomainNode{1},1), Coor(DomainNode{1},2), '.r');
hold on;
% plot(Coor(DomainNode{2},1), Coor(DomainNode{2},2), '.g');
% hold on;
% plot(Coor(DomainNode{3},1), Coor(DomainNode{3},2), '.b');
% hold on;
% plot(Coor(DomainNode{4},1), Coor(DomainNode{4},2), '.y');
% hold on;
axis equal;
%--------------------------------��֤C++���ɵ�d_node_pos�Ƿ���ȷ
d_dof = csvread('../../matrix/d_dof.csv');
d_node_pos = cell(4,1);
d_freenodesx = cell(4,1);
d_freenodesy = cell(4,1);
d_nodeid = cell(4,1);
figure(4);
patch(Domainx{1}',Domainy{1}','red','FaceAlpha',.3);
hold on;
patch(Domainx{2}',Domainy{2}','green','FaceAlpha',.3);
hold on;
patch(Domainx{3}',Domainy{3}','blue','FaceAlpha',.3);
hold on;
patch(Domainx{4}',Domainy{4}','yellow','FaceAlpha',.3);
hold on;
axis equal;
for i = 1:4
    d_freenodesx{i} = zeros(d_dof(i),1);
    d_freenodesy{i} = zeros(d_dof(i),1);
    d_nodeid{i} = zeros(d_dof(i),1);
    str = sprintf('../../matrix/d_node_pos[%d].csv', i-1);
    d_node_pos{i} =  csvread(str) + 1;
    str = sprintf('../../matrix/d_nodeid[%d].csv', i-1);
    d_nodeid{i} =  csvread(str) + 1;
    for n = 1:length(d_node_pos{i})
        if(d_node_pos{i}(n) <= d_dof(i))
            globalid = d_nodeid{i}(n);
            reorderid = d_node_pos{i}(n);
            d_freenodesx{i}(reorderid) = Coor(globalid,1);
            d_freenodesy{i}(reorderid) = Coor(globalid,2);
        end
    end
end
plot(d_freenodesx{1}, d_freenodesy{1}, '.r');
hold on;
plot(d_freenodesx{2}, d_freenodesy{2}, '.g');
hold on;
plot(d_freenodesx{3}, d_freenodesy{3}, '.b');
hold on;
plot(d_freenodesx{4}, d_freenodesy{4}, '.r');
hold on;
%--------------------------------����������Է���
At_real = csvread('../../matrix/At_real.csv');
figure(5);
plot(At_real,'ob');
hold on;
for i = 1:300
    str = sprintf('../../matrix/At_step%d.csv', i-1);
    At_step = csvread(str);
    str = sprintf('Step: %d', i);
    title(str);
    h = plot(At_step, '.r');
    pause(0.5);
    delete(h);
end
% % ------------------------------����0ֵ���ڵ�λ��
% At_real = csvread('../../matrix/At_real.csv');
% At_step = csvread('../../matrix/At_step99.csv');
% Atreal_zeronode = find(At_real==0);
% Atstep_zeronode = find(At_step==0);
% figure(5);
% plot(Coor(Atreal_zeronode,1), Coor(Atreal_zeronode,2), '.b');
% hold on;
% plot(Coor(Atstep_zeronode,1), Coor(Atstep_zeronode,2), '.r');