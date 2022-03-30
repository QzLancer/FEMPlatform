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
% %--------------------------------��ȡ��Ԫ����������
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
% %---------------------------------��ȡ�ڵ����������
% nDomain0Coor = Coor(nDomain==0,:);
% nDomain1Coor = Coor(nDomain==1,:);
% plot(nDomain0Coor(:,1),nDomain0Coor(:,2),'.b');
% hold on;
% axis equal;
% plot(nDomain1Coor(:,1),nDomain1Coor(:,2),'.r');
% hold on;
% %---------------------------------�ҳ������߽��ĵ�Ԫ
% %�������е�Ԫ���ҳ���Щ��Ԫ�����˲�ֹһ����Ľڵ�
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
% %��ȡ�����Ƴ���Щ��Ԫ
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
%�ҳ����ظ������Ľڵ�
% Domain0BoundNode = TriElement(Bound0Element,:);
% Domain1BoundNode = TriElement(Bound1Element,:);
% intersect(Domain0BoundNode,Domain1BoundNode);
%-----------------------------------�ҳ��߽�ڵ�
% %����ȫ����Ԫ��ȫ���ڵ㣬������ֽڵ��������뵥Ԫ������һ�£����Ϊ���紦�ڵ�
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