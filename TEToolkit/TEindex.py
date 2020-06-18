

"""Module Description

Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@author:  Ying Jin
@contact: yjin@cshl.edu
"""
import sys
import time
import logging
from math import ceil, floor
from TEToolkit.Constants import TEindex_BINSIZE

sys.setrecursionlimit(10000)


class Node:
    def __init__(self, start=-1, end=-1, name=-1, parent=None, left=None, right=None):
        self.__start = start
        self.__end = end
        self.__name = name  # idx in nameIDmap
        self.__namelist = {}
        self.left = left
        self.right = right
        self.balanceFactor = 0
        self.parent = parent
        self.isroot = False
        self.add(start, end, name)

    def add(self, start, end, name):
        if start in self.__namelist:
            self.__namelist[start].append((name, end))
        else:
            self.__namelist[start] = [(name, end)]

    def isRoot(self):
        return self.isroot

    def isLeftChild(self):
        return self.parent and self.parent.left == self

    def isRightChild(self):
        return  self.parent and self.parent.right == self

    def getStart(self):
        bin_startID = self.__start//TEindex_BINSIZE
        if self.__start == bin_startID * TEindex_BINSIZE:
            bin_startID -= 1
        return bin_startID

    def getEnd(self):
        bin_endID = self.__end//TEindex_BINSIZE
        return bin_endID

    def getName(self):
        return self.__name

    def overlaps(self, start, end):
        TEnamelist = []
        for s in sorted(self.__namelist.keys()):
            if s > end:
                break
            eles = self.__namelist[s]
            for name, e in eles:
                if start <= e and end >= s:
                    TEnamelist.append(name)

        return TEnamelist


class BinaryTree:
    def __init__(self):

        self.root = None
        self.size = 0

    def children_count(self):
        """
        Return the number of children

        @returns number of children: 0, 1, 2
        """
        cnt = 0
        if self.left:
            cnt += 1
        if self.right:
            cnt += 1
        return cnt

    def insert(self, start, end, name):
        if self.root:
            self.__insert(self.root, start, end, name)
        else:
            self.root = Node(start, end, name)
            self.root.isroot = True
        self.size = self.size + 1

    def __insert(self, node, start, end, name):
        """
        Insert new node with data

        @param data node data object to insert
        """
        root = node
        binstart = start//TEindex_BINSIZE
        if start == binstart * TEindex_BINSIZE:
            binstart -= 1

        if binstart == root.getStart():
            root.add(start, end, name)

        if binstart < root.getStart():
            if root.left is None:
                root.left = Node(start=start, end=end, name=name, parent=root)

                self.updateBalance(root.left)
            else:
                self.__insert(root.left, start, end, name)
        if binstart > root.getStart():
            if root.right is None:
                root.right = Node(start, end, name, parent=root)

                self.updateBalance(root.right)
            else:
                self.__insert(root.right, start, end, name)

    def updateBalance(self, node):
        if node.balanceFactor > 1 or node.balanceFactor < -1:
            self.rebalance(node)
            return

        if node.parent is not None:
            if node.isLeftChild():
                node.parent.balanceFactor += 1
            elif node.isRightChild():
                node.parent.balanceFactor -= 1

            if node.parent.balanceFactor != 0:
                self.updateBalance(node.parent)

    def rebalance(self, node):
        if node.balanceFactor < 0:
            if node.right.balanceFactor > 0:
                self.rotateRight(node.right)
                self.rotateLeft(node)
            else:
                self.rotateLeft(node)
        elif node.balanceFactor > 0:
            if node.left.balanceFactor < 0:
                self.rotateLeft(node.left)
                self.rotateRight(node)
            else :
                self.rotateRight(node)

    def rotateRight(self, oldRoot):
        newRoot = oldRoot.left
        oldRoot.left = newRoot.right
        if newRoot.right is not None:
            newRoot.right.parent = oldRoot
        newRoot.parent = oldRoot.parent

        if oldRoot.isRoot():
            self.root = newRoot
            newRoot.isroot = True
        else:
            if oldRoot.isLeftChild():
                oldRoot.parent.left = newRoot
            else:
                oldRoot.parent.right = newRoot

        newRoot.right = oldRoot
        oldRoot.parent = newRoot
        oldRoot.balanceFactor = oldRoot.balanceFactor - 1 - max(newRoot.balanceFactor, 0)
        newRoot.balanceFactor = newRoot.balanceFactor - 1 + min(oldRoot.balanceFactor, 0)

    def rotateLeft(self, oldRoot):
        newRoot = oldRoot.right
        oldRoot.right = newRoot.left
        if newRoot.left is not None:
            newRoot.left.parent = oldRoot
        newRoot.parent = oldRoot.parent

        if oldRoot.isRoot():
            self.root = newRoot
            newRoot.isroot = True
        else:
            if oldRoot.isLeftChild():
                oldRoot.parent.left = newRoot
            else:
                oldRoot.parent.right = newRoot

        newRoot.left = oldRoot
        oldRoot.parent = newRoot
        oldRoot.balanceFactor = oldRoot.balanceFactor + 1 - min(newRoot.balanceFactor, 0)
        newRoot.balanceFactor = newRoot.balanceFactor + 1 - max(oldRoot.balanceFactor, 0)

    # range query
    def lookup_r(self, start, end, node):

        if node is None:
            return None, None

        node_start = node.getStart()
        if end < node_start:
            return self.lookup_r(start, end, node.left)

        if start > node_start:
            return self.lookup_r(start, end, node.right)

        if start == node_start and end == node_start:
            return node, None

        if start == node_start and end > node_start:
            return node, self.lookup_p(end, node.right)

        if end == node_start and start < node_start:
            return self.lookup_p(start, node.left), node
        return None, None

    # point query using start point
    def lookup_p(self, start, node):
        """
        Lookup node containing data

        @param data node data object to look up
        @param parent node's parent
        @returns node and node's parent if found or None, None
        """
        if node is None:
            return None
        if start < node.getStart():
            if node.left is None:
                return None
            return self.lookup_p(start, node.left)

        elif start > node.getStart():
            if node.right is None:
                return None
            return self.lookup_p(start, node.right)
        else:
            return node


class TEfeatures:
    """index of TE annotations.
    """
    def __init__ (self):

        self.indexlist = {}
        self._length = []
        self._nameIDmap = []
        self._elements = []

    def getNames(self):
        names = []
        return self._nameIDmap

    def numInstances(self):
        return len(self._nameIDmap)

    def getElements(self):
        return self._elements

    def getStrand(self, idx):
        f_name = self._nameIDmap[idx]
        return f_name[len(f_name)-1]

    def getEleName(self, idx):
        full_name = None
        if idx >= len(self._nameIDmap) or idx < 0:
            return None
        else:
            full_name = self._nameIDmap[idx]
        if full_name is not None:
            pos = full_name.find(':')
            val = full_name[pos+1:(len(full_name)-2)]
            return val
        else:
            return None

    def getFullName(self, idx):
        if idx >= len(self._nameIDmap) or idx < 0:
            return None
        else:
            return self._nameIDmap[idx]

    def getLength(self, TE_name_idx):
        if TE_name_idx < len(self._length):
            return self._length[TE_name_idx]
        else:
            return -1

    def getFamilyID(self, chromosome, start, end):
        binID = start//TEindex_BINSIZE
        endbinID = end//TEindex_BINSIZE + 1

        if chromosome in self.indexlist:
            index = self.indexlist[chromosome]
            (node, RBnode) = index.lookup(binID, index.root, None, None)

            if node is not None and node.overlaps(binID, endbinID):
                full_name = (node.getName()).split(':')
                famid = full_name[2]
                return famid
            else:
                return None

        else:
            return None

    def findOvpTE(self, chrom, start, end):
        startbinID = start//TEindex_BINSIZE
        endbinID = end//TEindex_BINSIZE
        if start == startbinID * TEindex_BINSIZE:
            startbinID -= 1
        name_idx_list = []

        if chrom in self.indexlist:
            index = self.indexlist[chrom]
        else:
            return None

        (LBnode, RBnode) = index.lookup_r(startbinID, endbinID, index.root)

        if LBnode is not None:
            telist = LBnode.overlaps(start, end)
            name_idx_list.extend(telist)

        if RBnode is not None:
            telist = RBnode.overlaps(start, end)
            name_idx_list.extend(telist)

        return name_idx_list

    def TE_annotation(self, iv_seq):
        TEs = []
        for iv in iv_seq:
            chromosome = iv[0]
            start = iv[1]
            end = iv[2]
            strand = iv[3]
            name_idx_list = self.findOvpTE(chromosome, start, end)

            if name_idx_list is not None:
                for t in name_idx_list:
                    if strand != ".":  # stranded
                        if strand == self.getStrand(t):
                            if t not in TEs:
                                TEs.append(t)
                    else:  # not stranded
                        if t not in TEs:
                            TEs.append(t)

        return TEs

    # group by element
    def groupByEle(self, te_inst_counts):

        TEs = self.getElements()
        te_ele_counts = dict(list(zip(TEs, [0]*len(TEs))))

        for i in range(len(te_inst_counts)):
            ele_name = self.getEleName(i)

            if ele_name is None:
                sys.stderr.write("TE out of index boundary! \n")
                sys.exit(1)

            if ele_name in te_ele_counts:
                te_ele_counts[ele_name] += te_inst_counts[i]
            else:
                sys.stderr.write("TE inconsistency! "+ele_name+" \n")
                sys.exit(1)

        return te_ele_counts

    def build(self, filename):
        self.__srcfile = filename

        try:
            f = open(self.__srcfile, 'r')

        except:
            logging.error("cannot open such file %s ! \n" %(self.__srcfile))
            sys.exit(1)

        name_idx = 0
        linenum = 0
        for line in f:
            line = line.strip()

            if line.startswith("#"):
                continue

            items = line.split('\t')
            chrom = items[0]
            start = int(items[3])
            end = int(items[4])
            strand = items[6]
            items[8] = items[8].replace("; ", ";")
            desc = items[8].split(';')
            name = ""
            family_id = ""
            ele_id = ""
            class_id = ""
            tlen = end - start + 1
            linenum += 1

            for i in range(len(desc)):
                desc[i] = desc[i].replace("\"", "")
                pos = desc[i].find(" ")
                tid = desc[i][:pos]
                val = desc[i][pos+1:len(desc[i])]

                if tid == "gene_id":
                    ele_id = val
                if tid == "transcript_id":
                    name = val
                if tid == "family_id":
                    family_id = val
                if tid == "class_id":
                    class_id = val

            if ele_id == "" or name == "" or family_id == "" or class_id == "":
                sys.stderr.write(line+" \n")
                sys.stderr.write("TE GTF format error! There is no annotation at line %s. \n" % linenum)
                raise

            full_name = name+':'+ele_id+':'+family_id+':'+class_id+':'+strand
            ele_name = ele_id+':'+family_id+':'+class_id
            if ele_name not in self._elements:
                self._elements.append(ele_name)

            self._length.append(tlen)
            self._nameIDmap.append(full_name)

            if chrom in self.indexlist:
                index = self.indexlist[chrom]

                bin_startID = start//TEindex_BINSIZE
                bin_endID = end//TEindex_BINSIZE

                if start == bin_startID * TEindex_BINSIZE:
                    bin_startID -= 1
                while bin_startID <= bin_endID:
                    end_pos = min(end, (bin_startID+1) * TEindex_BINSIZE)
                    start_pos = max(start, bin_startID * TEindex_BINSIZE+1)

                    index.insert(start_pos, end_pos, name_idx)
                    bin_startID += 1

            else:
                index = BinaryTree()
                bin_startID = start//TEindex_BINSIZE
                bin_endID = end//TEindex_BINSIZE

                if start == bin_startID * TEindex_BINSIZE:
                    bin_startID -= 1

                while bin_startID <= bin_endID:
                    end_pos = min(end, (bin_startID+1) * TEindex_BINSIZE)
                    start_pos = max(start, bin_startID * TEindex_BINSIZE+1)
                    index.insert(start_pos, end_pos, name_idx)
                    bin_startID += 1

                self.indexlist[chrom] = index

            name_idx += 1

        f.close()


if __name__ == '__main__':
    try:

        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt! \n")
        sys.exit(0)
