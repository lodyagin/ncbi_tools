/*  parttree.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  parttree.c
*
* Author:  Vladimir Soussov
*   
* File Description:  taxonomy partial tree functions
*
*
* $Log: parttree.c,v $
* Revision 1.4  1998/07/23 18:24:30  soussov
* Added function for attaching subtrees
*
* Revision 1.3  1998/06/11 20:43:56  soussov
* kill some warnings
*
* Revision 1.2  1998/04/01 17:34:08  soussov
* changed tp include <>
*
* Revision 1.1  1998/02/10 20:11:57  soussov
* taxon2 related soft
*
 * Revision 1.1.1.1  1997/04/30  21:29:47  soussov
 * initial tree for taxon2
 *
*
*/

#include <txclient.h>

static VoidPtr ptree_getNode(TreeCursorPtr cursor, Uint4 format, Uint2Ptr size)
{
    Uint2 s;
    TXC_TreeNodePtr tnp;

    *size= 4;
    if(format == TREE_NODE_FMT_DISTANCE) return (VoidPtr)1;
    if(format == TREE_NODE_FMT_ICON) return (VoidPtr)-1;

    tnp= tree_getNodeData(cursor, &s);
    switch(format) {
    case TREE_NODE_FMT_LABEL: 
	*size= strlen(tnp->node_label);
	return tnp->node_label;
    case TREE_NODE_FMT_DEFAULT:
	*size= s;
	return tnp;
    default:
	return &tnp->tax_id;
    }
}
	

TreePtr tax_ptree_new(void)
{
    TreePtr ptree= tree_new();
    TreeCursorPtr cursor= tree_openCursor(ptree, NULL, NULL);
    TXC_TreeNodePtr tnp= MemNew(20);

    tnp->tax_id= 1;
    tnp->flags= 0;
    StringCpy(tnp->node_label, "root");

    tree_updateNodeData(cursor, tnp, 20);

    tree_setGetNodeFunc(ptree, ptree_getNode);

    MemFree(tnp);
    tree_closeCursor(cursor);
    return ptree;
}

/********************************************
 * Find id among children of current node and make this node current
 */
static Boolean known_id(TreeCursorPtr cursor, Int4 id)
{
    TXC_TreeNodePtr tnp;
    Uint2 s;

    if(tree_child(cursor)) {
	do {
	    tnp= tree_getNodeData(cursor, &s);
	    if(tnp->tax_id == id) return TRUE;
	}
	while(tree_sibling(cursor));
	tree_parent(cursor);
    }
    return FALSE;
}
		    
    
Boolean tax_ptree_addNode(TreePtr ptree, Int4 tax_id)
{
    Int4 lin_size, i;
    Uint2 s;
    TXC_TreeNodePtr* lineage;
    TreeCursorPtr cursor;
    TreeNodeId nid;

    lineage= txc_getLineage(tax_id, &lin_size);

    if((lineage == NULL) || (lin_size < 1)) return FALSE;

    cursor= tree_openCursor(ptree, NULL, NULL);
    if(cursor == NULL) return FALSE;
    
    for(i= lin_size; i--; ) {
	if(!known_id(cursor, lineage[i]->tax_id)) break;
	MemFree(lineage[i]);
    }
    
    if(i >= 0) {
    
	/* add nodes begining with this point */
	do {
	    s= 9 + StringLen(lineage[i]->node_label);
	    nid= tree_addChild(cursor, lineage[i], s);
	    tree_toNode(cursor, nid);
	    MemFree(lineage[i]);
	}
	while(i--);
    }

    tree_closeCursor(cursor);
    MemFree(lineage);
    return TRUE;
}

static Boolean ptree_toTaxId(TreeCursorPtr cursor, Int4 id)
{
    TXC_TreeNodePtr tnp;
    Uint2 s;

    tnp= tree_getNodeData(cursor, &s);
    if(tnp->tax_id == id) return TRUE;

    if(tree_child(cursor)) {
	do {
	    if(ptree_toTaxId(cursor, id)) return TRUE;
	}
	while(tree_sibling(cursor));

	tree_parent(cursor);
    }
    return FALSE;
}   

Boolean tax_ptree_toTaxId(TreeCursorPtr cursor, Int4 tax_id, Boolean search_in_subtree)
{
    TreeNodeId id= tree_getId(cursor);

    if(!search_in_subtree) {
	tree_root(cursor);
    }

    if(ptree_toTaxId(cursor, tax_id)) return TRUE;
    else {
	tree_toNode(cursor, id);
	return FALSE;
    }
}

Boolean tax_ptree_addChildren(TreeCursorPtr cursor)
{
    Int4 nof_children, i;
    Uint2 s;
    TXC_TreeNodePtr* child;
    TXC_TreeNodePtr my_node= tree_getNodeData(cursor, &s);
    TreeNodeId nid= tree_getId(cursor);
    Boolean need_remove= FALSE;

    if(my_node == NULL) return FALSE;

    child= txc_getChildren(my_node->tax_id, &nof_children);

    if((child == NULL) || (nof_children < 1)) return FALSE;

    if(tree_child(cursor)) {
	/* update all children which already exist in the tree */
	do {
	    my_node= tree_getNodeData(cursor, &s);
	    if(my_node != NULL) {
		for(i= 0; i < nof_children; i++) {
		    if((child[i] != NULL) && (child[i]->tax_id == my_node->tax_id)) {
			s= 9 + StringLen(child[i]->node_label);
			tree_updateNodeData(cursor, child[i], s);
			child[i]= MemFree(child[i]);
			break;
		    }
		}
		if(i >= nof_children) {
		    /* this node does not exists any more --> mark it for removing */
		    my_node->tax_id= 0;
		    need_remove= TRUE;
		}
	    }
	}
	while(tree_sibling(cursor));

	tree_parent(cursor);

	while(need_remove) {
	    /* remove all marked nodes */
	    need_remove= FALSE;
	    if(tree_child(cursor)) {
	    
		do {
		    my_node= tree_getNodeData(cursor, &s);
		    if((my_node != NULL) && (my_node->tax_id == 0)) {
			need_remove= TRUE;
			tree_delSubtree(cursor);
			break;
		    }
		}
		while (tree_sibling(cursor));
	    }
	}
    }

    /* add new children */
    tree_toNode(cursor, nid);
		    
    for(i= 0; i < nof_children; i++) {
	if(child[i] != NULL) {
	    s= 9 + StringLen(child[i]->node_label);
	    tree_addChild(cursor, child[i], s);
	    child[i]= MemFree(child[i]);
	}
    }
	
    MemFree(child);
    return TRUE;
}

Boolean tax_ptree_addSubtree(TreeCursorPtr cursor)
{
    Int4 nof_children, i;
    Uint2 s;
    TXC_TreeNodePtr* child;
    TXC_TreeNodePtr my_node= tree_getNodeData(cursor, &s);
    TreeNodeId nid= tree_getId(cursor);
    Boolean need_remove= FALSE;

    if(my_node == NULL) return FALSE;

    child= txc_getChildren(-(my_node->tax_id), &nof_children);

    if((child == NULL) || (nof_children < 1)) return FALSE;

    while(tree_child(cursor)) {
	tree_delSubtree(cursor);
    }

    /* add new children */
    tree_toNode(cursor, nid);
		    
    for(i= 0; i < nof_children; i++) {
	if(child[i] != NULL) {
	    if((child[i]->flags & TXC_SUFFIX) != 0) {
		tree_root(cursor);
		ptree_toTaxId(cursor, child[i]->tax_id);
	    }
	    else {
		s= 9 + StringLen(child[i]->node_label);
		tree_addChild(cursor, child[i], s);
	    }
	    child[i]= MemFree(child[i]);
	}
    }
	
    MemFree(child);
    return TRUE;
}
