% metis只能实现非重叠子域区域分解，通过metis输出的数据格式，得到重叠子域的区域分解
% 4个分区
% metis分网提供两组信息：1、单元所在非重叠子域编号 2、节点所在子域
% 1、保存每个节点周边的单元编号
% 2、遍历节点所在的子域，将节点周围的单元都标记该节点所在子域编号，得到单元所在重叠子域编号
% 3、遍历子域中的所有单元，如果单元的节点所在子域和单元所在子域不一样，为边界节点
% 遍历全部节点，检索包含该节点的全部单元，即可得到重叠版本的区域分解
% By QzLancer
% 2022/3/28

%---------------------------------读取文件
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('model3656.mphtxt');
[fileID] = fopen('ContactLinear_7172.metis.epart.4');
eDomain = fscanf(fileID,'%d\n');
[fileID] = fopen('ContactLinear_7172.metis.npart.4');
nDomain = fscanf(fileID,'%d\n');
%--------------------------------提取节点周围的全部单元
EleSizeinNode = zeros(length(Coor),1);
GlobalEleIDinNode = zeros(length(Coor), 10);
for e = 1 : length(TriElement)
    for n = 1 :3
        nodeID = TriElement(e,n);
        EleSizeinNode(nodeID) = EleSizeinNode(nodeID) + 1;
        GlobalEleIDinNode(nodeID, EleSizeinNode(nodeID)) = e;
    end
end
%--------------------------------遍历节点所在的子域，将相关单元保存epartTable中
ePartTable = zeros(length(TriElement), 4);
for n = 1 : length(Coor)
    domain = nDomain(n) + 1;
    for e = 1 : EleSizeinNode(n)
        GlobalEleID = GlobalEleIDinNode(n, e);
        ePartTable(GlobalEleID, domain) = 1;
    end
end
%----------------------------------读取单元分区并绘制
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
%----------------------------------检索每个单元的边界节点，保存在boundaryTable中
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
%----------------------------------绘制boundaryTable
DomainNode = cell(4,1); 
for i = 1:4
   DomainNode{i} = find(boundaryTable(:,i)); 
end
plot(Coor(DomainNode{1},1), Coor(DomainNode{1},2), '.r');
plot(Coor(DomainNode{2},1), Coor(DomainNode{2},2), '.g');
plot(Coor(DomainNode{3},1), Coor(DomainNode{3},2), '.b');
plot(Coor(DomainNode{4},1), Coor(DomainNode{4},2), '.y');
%----------------------------------验证C++生成的eparttable和boundaryTable是否正确
figure(2);
ePartTable1 = csvread('../../matrix/eparttable.csv');
ePartTable1 = ePartTable1(:,1:length(ePartTable))';
boundaryTable1 = csvread('../../matrix/boundarytable.csv');
boundaryTable1 = boundaryTable1(:,1:length(boundaryTable))';
DomainElement = cell(4,1);
Domainx = cell(4,1);
Domainy = cell(4,1);
for i = 1:4
    DomainElement{i} = find(ePartTable1(:, i));
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
DomainNode = cell(4,1); 
for i = 1:4
   DomainNode{i} = find(boundaryTable1(:,i)); 
end
plot(Coor(DomainNode{1},1), Coor(DomainNode{1},2), '.r');
plot(Coor(DomainNode{2},1), Coor(DomainNode{2},2), '.g');
plot(Coor(DomainNode{3},1), Coor(DomainNode{3},2), '.b');
plot(Coor(DomainNode{4},1), Coor(DomainNode{4},2), '.y');
% %--------------------------------读取单元分区并绘制
% Domain0Element = TriElement(eDomain==0,:);
% Domain1Element = TriElement(eDomain==1,:);
% Domain2Element = TriElement(eDomain==2,:);
% Domain3Element = TriElement(eDomain==3,:);
% Domain0x(:,1) = Coor(Domain0Element(:,1),1);
% Domain0x(:,2) = Coor(Domain0Element(:,2),1);
% Domain0x(:,3) = Coor(Domain0Element(:,3),1);
% Domain0y(:,1) = Coor(Domain0Element(:,1),2);
% Domain0y(:,2) = Coor(Domain0Element(:,2),2);
% Domain0y(:,3) = Coor(Domain0Element(:,3),2);
% patch(Domain0x',Domain0y','red','FaceAlpha',.3);
% hold on;
% Domain1x(:,1) = Coor(Domain1Element(:,1),1);
% Domain1x(:,2) = Coor(Domain1Element(:,2),1);
% Domain1x(:,3) = Coor(Domain1Element(:,3),1);
% Domain1y(:,1) = Coor(Domain1Element(:,1),2);
% Domain1y(:,2) = Coor(Domain1Element(:,2),2);
% Domain1y(:,3) = Coor(Domain1Element(:,3),2);
% patch(Domain1x',Domain1y','blue','FaceAlpha',.3);
% hold on;
% Domain2x(:,1) = Coor(Domain2Element(:,1),1);
% Domain2x(:,2) = Coor(Domain2Element(:,2),1);
% Domain2x(:,3) = Coor(Domain2Element(:,3),1);
% Domain2y(:,1) = Coor(Domain2Element(:,1),2);
% Domain2y(:,2) = Coor(Domain2Element(:,2),2);
% Domain2y(:,3) = Coor(Domain2Element(:,3),2);
% patch(Domain2x',Domain2y','yellow','FaceAlpha',.3);
% hold on;
% Domain3x(:,1) = Coor(Domain3Element(:,1),1);
% Domain3x(:,2) = Coor(Domain3Element(:,2),1);
% Domain3x(:,3) = Coor(Domain3Element(:,3),1);
% Domain3y(:,1) = Coor(Domain3Element(:,1),2);
% Domain3y(:,2) = Coor(Domain3Element(:,2),2);
% Domain3y(:,3) = Coor(Domain3Element(:,3),2);
% patch(Domain3x',Domain3y','green','FaceAlpha',.3);
% hold on;
% %---------------------------------读取节点分区并绘制
% nDomain0Coor = Coor(nDomain==0,:);
% nDomain1Coor = Coor(nDomain==1,:);
% plot(nDomain0Coor(:,1),nDomain0Coor(:,2),'.b');
% hold on;
% axis equal;
% plot(nDomain1Coor(:,1),nDomain1Coor(:,2),'.r');
% hold on;
% %---------------------------------找出包含边界点的单元
% %遍历所有单元，找出那些单元包含了不止一个域的节点
% eleDomainNode = nDomain(TriElement);
% j = 1;
% k = 1
% for i = 1:length(TriElement)
%     if length(unique(eleDomainNode(i,:))) ~= 1
%         if eDomain(i) == 0
%             Bound0Element(j,1) = i;
%             j = j+1;
%         else
%             Bound1Element(k,1) = i;
%             k = k+1;
%         end
%     end
% end
% %提取并绘制出这些单元
% DomainB0Node = TriElement(Bound0Element,:);
% eDomainB0x(:,1) = Coor(DomainB0Node(:,1),1);
% eDomainB0x(:,2) = Coor(DomainB0Node(:,2),1);
% eDomainB0x(:,3) = Coor(DomainB0Node(:,3),1);
% eDomainB0y(:,1) = Coor(DomainB0Node(:,1),2);
% eDomainB0y(:,2) = Coor(DomainB0Node(:,2),2);
% eDomainB0y(:,3) = Coor(DomainB0Node(:,3),2);
% patch(eDomainB0x',eDomainB0y','blue','FaceAlpha',.3);
% hold on;
% DomainB1Node = TriElement(Bound1Element,:);
% eDomainB1x(:,1) = Coor(DomainB1Node(:,1),1);
% eDomainB1x(:,2) = Coor(DomainB1Node(:,2),1);
% eDomainB1x(:,3) = Coor(DomainB1Node(:,3),1);
% eDomainB1y(:,1) = Coor(DomainB1Node(:,1),2);
% eDomainB1y(:,2) = Coor(DomainB1Node(:,2),2);
% eDomainB1y(:,3) = Coor(DomainB1Node(:,3),2);
% patch(eDomainB1x',eDomainB1y','red','FaceAlpha',.3);
% hold on;
%找出被重复遍历的节点
% Domain0BoundNode = TriElement(Bound0Element,:);
% Domain1BoundNode = TriElement(Bound1Element,:);
% intersect(Domain0BoundNode,Domain1BoundNode);
%-----------------------------------找出边界节点
% %遍历全部单元的全部节点，如果出现节点所在域与单元所在域不一致，标记为交界处节点
% eleDomainNode = nDomain(TriElement);
% k = 1; 
% for i = 1:length(eDomain)
%     for j = 1:3
%         if eleDomainNode(i,j) ~= eDomain(i)
%             DomainBNode(k,1) = TriElement(i,j);
%             k = k+1;
%         end
%     end
% end
% plot(Coor(DomainBNode,1),Coor(DomainBNode,2),'.g');