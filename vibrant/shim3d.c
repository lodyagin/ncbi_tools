/*   shim3d.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  shim3d.c
*
* Author:  Lewis Geer
*
* Version Creation Date:   1/26/99
*
* $Revision: 6.1 $
*
* File Description: Shim functions to replace Viewer3D with OpenGL
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: shim3d.c,v $
* Revision 6.1  1999/04/06 14:23:28  lewisg
* add opengl replacement for viewer3d
*
*
*/


#ifdef _OPENGL

/*
*  Include the GL dependencies.  GL has its own typedef's for basic types, just like the toolkit.
*  If you get warnings about type mismatch, this should be investigated.  The GL typedef's can't
*  be included in general toolkit code because of the windows.h dependency for WIN32 which
*  causes all sorts of name collisions.
*/

#ifdef WIN32  /* braindead windows dependency */
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>

/* 
*  the following 2 includes are a subset of vibrant.h, since vibdefns.h
*  and vibforms.h confict with windows.h 
*/

#include <vibtypes.h>
#include <vibprocs.h>

#include <math.h>
#include <shim3d.h>


#define OGL_SQR(oglx) ((oglx)*(oglx))

typedef struct _TOGL_Layers
/* this struct contains the information used to manage the different layers of the display */
{
    GLuint FirstLayer;
    GLuint LastLayer;
    GLuint SelectedLayer;
    Nlm_Boolean IsOn[OGLMAXLAYERS];
} TOGL_Layers;


/*
*   Various helper functions used in drawing 
*/

FloatHi * OGL_CrossProduct( Nlm_FloatHi * v1, Nlm_FloatHi * v2)
{
    Nlm_FloatHi * RetVal;
    
    if( v1 == NULL || v2 == NULL) return;
    RetVal = MemNew(3 * sizeof(Nlm_FloatHi));
    if (RetVal == NULL) return NULL;
    
    RetVal[0] = v1[1]*v2[2] - v1[2]*v2[1];
    RetVal[1] = v1[2]*v2[0] - v1[0]*v2[2];
    RetVal[2] = v1[0]*v2[1] - v1[1]*v2[0];
    
    return RetVal;
}


void OGL_Normalize( Nlm_FloatHi * v)
/* normalize a vector */
{
    Nlm_FloatHi Length;
    
    if (v == NULL) return;
    Length = sqrt(OGL_SQR(v[0])+OGL_SQR(v[1])+OGL_SQR(v[2]));
    v[0] /= Length;
    v[1] /= Length;
    v[2] /= Length;
}


FloatHi * OGL_MakeNormal(Nlm_FloatHi * origin, Nlm_FloatHi * v1, Nlm_FloatHi * v2)
/* creates a normal to the given 3 vertices */
{
    Nlm_FloatHi Vector1[3], Vector2[3], * RetValue;
    
    if(origin == NULL || v1 == NULL || v2 == NULL) return NULL;
    Vector1[0] = v1[0] - origin[0];
    Vector1[1] = v1[1] - origin[1];
    Vector1[2] = v1[2] - origin[2];
    
    Vector2[0] = v2[0] - origin[0];
    Vector2[1] = v2[1] - origin[1];
    Vector2[2] = v2[2] - origin[2];
    
    RetValue = OGL_CrossProduct(Vector1, Vector2);
    OGL_Normalize(RetValue);
    return RetValue;
}


static void ColorCell2Array(GLfloat * array, TOGL_ColorCell * color)
/* copies a color cell to a GL array */
{
    if (array == NULL || color == NULL) return;
    array[0] = (GLfloat)(color->rgb[0]/255.0);
    array[1] = (GLfloat)(color->rgb[1]/255.0);
    array[2] = (GLfloat)(color->rgb[2]/255.0);
}
    

/*
 *  Functions used to draw various primitives
 */

void OGL_AddQuad3D(TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                   Nlm_FloatHi * v1, Nlm_FloatHi * v2, Nlm_FloatHi * v3, Nlm_FloatHi * v4)
/* draws a quadralateral with the 4 given vertices of form double v1[3] */
{
    Nlm_FloatHi * Normal;
    GLfloat Ambient[] = { 0.0, 0.0, 0.0, 1.0 };  /* for RGBA mode */
    GLfloat Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };  /* for RGBA mode */
    ValNodePtr PaletteIndex;
    GLfloat ColorMap[3];  /* index color map */ 
    
    
    if (v1 == NULL || v2 == NULL || v3 == NULL || v4 == NULL || OGL_Data == NULL ||
        color == NULL) return;
    
    if(OGL_Data->IndexMode) {
        PaletteIndex = OGL_SearchPaletteIndex (OGL_Data->PaletteIndex, color);
        if(!PaletteIndex) return;
        ColorMap[0] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->Begin;
        ColorMap[1] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[2] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        glMaterialfv(GL_FRONT_AND_BACK, GL_COLOR_INDEXES, ColorMap);
    }
    else {
        ColorCell2Array(Diffuse, color);
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, Ambient);
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, Diffuse);
    }
    
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    
    glBegin(GL_QUADS);
    Normal = OGL_MakeNormal(v1, v2, v4);
    if(Normal != NULL) glNormal3dv(Normal);
    glVertex3dv(v1);
    glVertex3dv(v2);
    glVertex3dv(v3);
    glVertex3dv(v4);
    glEnd();

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}


void OGL_AddCylinder3D (TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                        Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                        Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2, Nlm_FloatHi radius)
                        /* create a cylinder with given endcaps and radius */
{
    GLfloat Ambient[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    ValNodePtr PaletteIndex;
    GLfloat ColorMap[3];  /* index color map */ 
    GLdouble Rotate[16];
    GLdouble Translate[3];
    Nlm_Int4 iCount;
    GLdouble length;
    GLUquadricObj * qobj;
    GLdouble a, b, c;  /* the normal z */
    GLdouble yy2, yy3;  /* the normal y */
    GLdouble xx1, xx2, xx3;  /* the normal x */
    Nlm_Int4 slices;
    
    if(OGL_Data == NULL || color == NULL) return;
    for (iCount=0; iCount<16; iCount++) Rotate[iCount]= 0.0;
    Rotate[15] = 1;  /* identity */
    
    Translate[0] = (x1+x2)/2;
    Translate[1] = (y1+y2)/2;
    Translate[2] = (z1+z2)/2;
    
    length = sqrt(OGL_SQR(x1-x2)+OGL_SQR(y1-y2)+OGL_SQR(z1-z2))/2.0;
    
    /* create the normal z */
    a = (x1-Translate[0])/length;
    b = (y1-Translate[1])/length;
    c = (z1-Translate[2])/length;
    
    /* create the normal y */
    
    yy2 = sqrt(1.0/(1.0+OGL_SQR(b)/OGL_SQR(c)));
    yy3 = -(b/c)*yy2;
    
    /* create the normal x */
    
    xx2 = sqrt(1.0/(pow(c,4.0)/(OGL_SQR(a)*OGL_SQR(b))+2.0*OGL_SQR(c)/OGL_SQR(a)
        + OGL_SQR(b)/OGL_SQR(a) + 1 + OGL_SQR(c)/OGL_SQR(b)));
    xx3 = xx2*c/b;
    xx1 = (-OGL_SQR(c)/(a*b) - b/a) * xx2;
    
    
    /* now use the normals to make the rotation matrix */
    
    Rotate[0] = xx1;
    Rotate[1] = xx2;
    Rotate[2] = xx3;
    
    Rotate[4] = 0.0;
    Rotate[5] = yy2;
    Rotate[6] = yy3;
    
    Rotate[8] = a;
    Rotate[9] = b;
    Rotate[10] = c;
        
    if(OGL_Data->IndexMode) {
        PaletteIndex = OGL_SearchPaletteIndex (OGL_Data->PaletteIndex, color);
        if(!PaletteIndex) return;
        ColorMap[0] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->Begin;
        ColorMap[1] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[2] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        glMaterialfv(GL_FRONT, GL_COLOR_INDEXES, ColorMap);
    }
    else {
        ColorCell2Array(Diffuse, color);
        glMaterialfv( GL_FRONT, GL_AMBIENT, Ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, Diffuse);
    }

    glPushMatrix();
    /* the rotation/translation duplicates work in algorend.c */
    glTranslated(Translate[0], Translate[1], Translate[2]);
    glMultMatrixd(Rotate);
    glTranslated(0.0, 0.0, -length);
    qobj = gluNewQuadric();
    if(qobj == NULL) return;
    gluQuadricDrawStyle(qobj, GLU_FILL);
    gluQuadricNormals(qobj, GLU_SMOOTH);
    slices = (Int4) (radius*6);
    if (slices < 6) slices = 6;
    gluCylinder(qobj, radius, radius, length*2.0 , slices, 1);
    glPopMatrix();
    gluDeleteQuadric(qobj);
}


void OGL_AddLine3D (TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                    Nlm_FloatHi x1, Nlm_FloatHi y1, Nlm_FloatHi z1,
                    Nlm_FloatHi x2, Nlm_FloatHi y2, Nlm_FloatHi z2)
                    /* draw a single line */
{
    GLfloat Ambient[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    ValNodePtr PaletteIndex;
    GLfloat ColorMap[3];  /* index color map */   
    
    if(OGL_Data == NULL || color == NULL) return;
    if(OGL_Data->IndexMode) {
        PaletteIndex = OGL_SearchPaletteIndex (OGL_Data->PaletteIndex, color);
        if(!PaletteIndex) return;
        ColorMap[0] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[1] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[2] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        glMaterialfv(GL_FRONT_AND_BACK, GL_COLOR_INDEXES, ColorMap);
    }
    else {
        ColorCell2Array(Diffuse, color);
        ColorCell2Array(Ambient, color);
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, Ambient);
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, Diffuse);
    }

    glBegin(GL_LINES);
    glVertex3d(x1, y1, z1);
    glVertex3d(x2, y2, z2);
    glEnd();     
}


void OGL_AddSphere3D (TOGL_Data * OGL_Data, TOGL_ColorCell * color,
                      Nlm_FloatHi x, Nlm_FloatHi y, Nlm_FloatHi z, Nlm_FloatHi radius)
                      /* draws a sphere */
{
    GLUquadricObj * qobj;   
    Nlm_Int4 slices, stacks;
    ValNodePtr PaletteIndex;
    GLfloat ColorMap[3];  /* index color map */ 
    GLfloat Ambient[] = { 0.0, 0.0, 0.0, 1.0 };  /* for RGBA mode */
    GLfloat Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };  /* for RGBA mode */
    
    if(OGL_Data == NULL || color == NULL) return;
    if(OGL_Data->IndexMode) {
        PaletteIndex = OGL_SearchPaletteIndex (OGL_Data->PaletteIndex, color);
        if(!PaletteIndex) return;
        ColorMap[0] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->Begin;
        ColorMap[1] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[2] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        glMaterialfv(GL_FRONT, GL_COLOR_INDEXES, ColorMap);
    }
    else {
        ColorCell2Array(Diffuse, color);
        glMaterialfv( GL_FRONT, GL_AMBIENT, Ambient);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, Diffuse);
    }


    glPushMatrix();
    glTranslated(x, y, z);
    qobj = gluNewQuadric();
    if (qobj == NULL) return;
    gluQuadricDrawStyle(qobj, GLU_FILL);
    gluQuadricNormals(qobj, GLU_SMOOTH);
    slices = (Int4)(radius*6);
    if (slices < 8) slices = 8;
    stacks = (Int4)(radius*3);
    if(stacks < 4) stacks = 4;
    gluSphere(qobj, radius, slices, stacks);
    glPopMatrix();
    gluDeleteQuadric(qobj);
}


void OGL_AddText3D (TOGL_Data * OGL_Data, TOGL_ColorCell * color, Nlm_CharPtr string,
                    Nlm_FloatHi x, Nlm_FloatHi y, Nlm_FloatHi z, Nlm_Int2 flags)
{
    ValNodePtr PaletteIndex;
    GLfloat ColorMap[3];  /* index color map */ 
    Nlm_Int4 Length, i;
    GLfloat Ambient[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };

    if(OGL_Data == NULL || color == NULL || string == NULL) return;
    if(OGL_Data->IndexMode) {
        PaletteIndex = OGL_SearchPaletteIndex (OGL_Data->PaletteIndex, color);
        if(!PaletteIndex) return;
        ColorMap[0] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[1] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        ColorMap[2] = (GLfloat)((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End;
        glMaterialfv(GL_FRONT_AND_BACK, GL_COLOR_INDEXES, ColorMap);
    }
    else {
        ColorCell2Array(Diffuse, color);
        ColorCell2Array(Ambient, color);
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, Ambient);
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, Diffuse);
    }

    
    glListBase (OGLFONTBASE); 
    
    Length = Nlm_StrLen(string);
    glRasterPos3d(x, y, z);
    if(flags & OGLTEXT3D_CENTER) {
        for ( i = 0; i < Length/2; i++) glBitmap((GLsizei)OGL_Data->SpaceWidth, (GLsizei)OGL_Data->SpaceHeight,
            0.0, 0.0, -(GLfloat)OGL_Data->SpaceWidth, 0.0, OGL_Data->Space);
    }
    if(flags & OGLTEXT3D_MIDDLE) glBitmap((GLsizei)OGL_Data->SpaceWidth, (GLsizei)OGL_Data->SpaceHeight,
        0.0f, 0.0f, 0.0f, (GLfloat)(-OGL_Data->SpaceHeight/2.0), OGL_Data->Space);
    
    glCallLists (Length, GL_UNSIGNED_BYTE, string);
    
    glListBase (0);
}


/*
*   Functions used to manage display lists
*/

void OGL_Start(TOGL_Data * OGL_Data, Nlm_Int1 List)
/* begin a display list */
{
    if (OGL_Data == NULL) return;
    if (List >= OGLMAXLAYERS) return;
    glNewList(List + OGL_Data->Layers->FirstLayer, GL_COMPILE);
    OGL_SetLayer(OGL_Data, List + OGL_Data->Layers->FirstLayer, TRUE);
}

void OGL_End()
/* end a display list */
{
    glEndList();
}


void OGL_SetLayers(TOGL_Data * OGL_Data, Nlm_Boolean Status)
/* set the status of all the layers */
{
    Nlm_Int4 i;
    
    if (OGL_Data == NULL) return;
    for (i = 0; i < OGLMAXLAYERS; i++) 
        OGL_Data->Layers->IsOn[i] = Status;
    return;
}

void OGL_SetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i, Nlm_Boolean Status)
/* set the status of a particular layer -- is it on or off? */
{
    if(OGL_Data == NULL) return;
    OGL_Data->Layers->IsOn[i - OGL_Data->Layers->FirstLayer] = Status;
    return;
}

Nlm_Boolean OGL_GetLayer(TOGL_Data * OGL_Data, Nlm_Int4 i)
/* return layer status */
{
    if(OGL_Data == NULL) return;
    return OGL_Data->Layers->IsOn[i - OGL_Data->Layers->FirstLayer];
}


void OGL_SetLayerTop3D(TOGL_Data * OGL_Data, Nlm_Int4 TopLayer)
/* set the highest value layer used */
{
    if(OGL_Data == NULL) return;
    OGL_Data->Layers->LastLayer = TopLayer + OGL_Data->Layers->FirstLayer;
    return;
}

void OGL_AllLayerOnProc(TOGL_Data * OGL_Data)
/* turn on all used layers */
{
    if(OGL_Data == NULL) return;
    OGL_Data->Layers->SelectedLayer = 0;
    return;
}

void OGL_RewindLayerProc(TOGL_Data * OGL_Data)
/* rewind to the first layer */
{
    if(OGL_Data == NULL) return;
    OGL_Data->Layers->SelectedLayer = OGL_Data->Layers->FirstLayer;
    return;
}


void OGL_PrevLayerProc(TOGL_Data * OGL_Data)
/* go back to the previous layer */
{
    if(OGL_Data == NULL) return;
    if(OGL_Data->Layers->SelectedLayer) {
        if(OGL_Data->Layers->SelectedLayer > OGL_Data->Layers->FirstLayer)
            OGL_Data->Layers->SelectedLayer--;
    }
    else OGL_Data->Layers->SelectedLayer = OGL_Data->Layers->FirstLayer;
    
    return;
}


void OGL_NextLayerProc(TOGL_Data * OGL_Data)
/* go to the next layer */
{
    if(OGL_Data == NULL) return;
    if(OGL_Data->Layers->SelectedLayer) {
        if(OGL_Data->Layers->SelectedLayer < OGL_Data->Layers->LastLayer)
            OGL_Data->Layers->SelectedLayer++;
    }
    else OGL_Data->Layers->SelectedLayer = OGL_Data->Layers->FirstLayer;
    
    return;
}


void OGL_Play(TOGL_Data * OGL_Data)
/* used to flip through layers in endless loop */
{
    if(OGL_Data == NULL) return;
    if(OGL_Data->Layers->SelectedLayer) {
        if(OGL_Data->Layers->SelectedLayer < OGL_Data->Layers->LastLayer)
            OGL_Data->Layers->SelectedLayer++;
        else OGL_Data->Layers->SelectedLayer = OGL_Data->Layers->FirstLayer;
    }
    else OGL_Data->Layers->SelectedLayer = OGL_Data->Layers->FirstLayer;
    
    return;
}


/*
 *  Color manipulation functions
 */


ValNodePtr OGL_SearchPaletteIndex (ValNodePtr PaletteIndex, TOGL_ColorCell * ColorCell)
{
    if(PaletteIndex == NULL || ColorCell == NULL) return NULL;
    for(; PaletteIndex; PaletteIndex = PaletteIndex->next)
        if(((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->ColorCell.rgb[0] == ColorCell->rgb[0] &&
            ((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->ColorCell.rgb[1] == ColorCell->rgb[1] &&
            ((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->ColorCell.rgb[2] == ColorCell->rgb[2])
        {
            return PaletteIndex;
        }
        return NULL;
}


/*
 *  The mouse and rotation/translation/zoom code
 */

static TOGL_Data * MAToOGL(MAPtr ma)
/* extracts OGL_Data out of the extra pointer in the mouse data */
{
    return (TOGL_Data *)MA_GetExtra( ma );
}


static TOGL_Data * PanelToOGL(Nlm_PaneL panel)
/* extract OGL_Data out of the panel data */
{
    MAPtr ma;
    Nlm_GetPanelExtra(panel, &ma);
    
    return MAToOGL( ma );
}


static void OGL_DrawViewer3D_CB(Nlm_PaneL panel)
/* callback */
{
    OGL_DrawViewer3D(PanelToOGL( panel ) );
}


/*
 *  Move
 */

static void OGL_Move3D(TOGL_Data * OGL_Data, Nlm_Int2 dx, Nlm_Int2 dy)
/* this moves the structure by manipulation the viewpoint used in gluLookAt() */
{
    Nlm_FloatLo HeightMaxSize, WidthMaxSize;
 
    if(OGL_Data == NULL) return;
    if(OGL_Data->Width>OGL_Data->Height) {
        HeightMaxSize = (Nlm_FloatLo)(1.5*OGL_Data->MaxSize);
        WidthMaxSize = HeightMaxSize*(OGL_Data->Width/OGL_Data->Height);
    }
    else {
        WidthMaxSize = (Nlm_FloatLo)(1.5*OGL_Data->MaxSize);
        HeightMaxSize = WidthMaxSize*(OGL_Data->Width/OGL_Data->Height);
    }
    
    OGL_Data->ViewPoint[0] -= (FloatLo)dx/OGL_Data->Width*WidthMaxSize;
    OGL_Data->ViewPoint[1] += (FloatLo)dy/OGL_Data->Height*HeightMaxSize;
    
}


/*
 *  Zoom
 */


static void OGL_Zoom3D(TOGL_Data * OGL_Data, Nlm_Int2 x1, Nlm_Int2 y1, Nlm_Int2 x2, Nlm_Int2 y2)
{
    Nlm_FloatLo Zoom;
    
    if(OGL_Data == NULL) return;
    OGL_Move3D(OGL_Data, (Nlm_Int2)(OGL_Data->Width/2.0 - (x2+x1)/2.0),
        (Nlm_Int2)(OGL_Data->Height/2.0 - (y2+y1)/2.0));
    
    if(abs(x2-x1) > abs(y2-y1)) Zoom = abs(y2-y1)/(FloatLo)OGL_Data->Height;
    else Zoom = abs(x2-x1)/(FloatLo)OGL_Data->Width;
    
    OGL_Data->ViewPoint[3] *= Zoom;
    
}


/*
 *  Rotation
 */


typedef enum
{
    ROTATE_X,
        ROTATE_Y,
        ROTATE_Z
}
OGL_enumRotate3D;

typedef struct
{
    OGL_enumRotate3D H; /* horizontal dragging */
    OGL_enumRotate3D V; /* vertical   dragging */
}
OGL_RotatePivots3D, * OGL_RotatePivots3DPtr;


static void OGL_Rotate(TOGL_Data * OGL_Data, Nlm_Int4 dAngle, OGL_enumRotate3D pivot)
{
    
    if ( !dAngle )
        return; 
    if(OGL_Data == NULL) return;

    glLoadIdentity();
    
    switch ( pivot )
    {
    case ROTATE_X:
        glRotatef((GLfloat)(dAngle), 1.0f, 0.0f, 0.0f);
        break;
    case ROTATE_Y:
        glRotatef((GLfloat)(dAngle), 0.0f, 1.0f, 0.0f);
        break;
    case ROTATE_Z:
        glRotatef((GLfloat)(dAngle), 0.0f, 0.0f, 1.0f);
        break;
    }
    
    glMultMatrixd((GLdouble *)OGL_Data->ModelMatrix);
    
    glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *)OGL_Data->ModelMatrix); 
    
}



static void OGL_ViewerRotate(Nlm_SlatE panel, Nlm_Int4 delta,
                             OGL_enumRotate3D pivot,
                             Nlm_Boolean adjust_scrollbar, Nlm_Boolean redraw)
{
    TOGL_Data * OGL_Data;

    if (panel == NULL  ||
        !Nlm_Visible(panel)  ||  !Nlm_AllParentsVisible(panel))
        return;
    
    OGL_Data = PanelToOGL( (Nlm_PaneL)panel );
    
    if ( adjust_scrollbar )
    { /* adjust the relevant rotation scrollbar, if any */
        Nlm_BaR sbar = NULL;
        switch ( pivot )
        {
        case ROTATE_X:
            sbar = Nlm_GetSlateVScrollBar( panel );
            break;
        case ROTATE_Y:
            sbar = Nlm_GetSlateHScrollBar( panel );
            break;
        case ROTATE_Z:
            sbar = OGL_Data->Z_rotate;
            break;
        }
        
        if ( sbar )
        {
            Nlm_ResetClip();
            Nlm_CorrectBarValue(sbar, (Nlm_GetValue(sbar) + delta + 360) % 360);
        }
    }
    
    /* the coordination transformation (rotation) */
    if (pivot == ROTATE_X)
        delta = -delta;
    OGL_Rotate(OGL_Data, delta, pivot);
    
    if ( redraw )
    { /* Draw the viewer */
        OGL_DrawViewer3D(OGL_Data);
    }
    
    /*  Nlm_DiagReset();*/
}



/* scrollbar callbacks */
static void OGL_ViewerVScrollProc(Nlm_BaR sb, Nlm_SlatE viewer, Nlm_Int2 newval, Nlm_Int2 oldval)
{
    OGL_ViewerRotate(viewer, newval - oldval, ROTATE_X, FALSE, TRUE);
}

static void OGL_ViewerHScrollProc(Nlm_BaR sb, Nlm_SlatE viewer, Nlm_Int2 newval, Nlm_Int2 oldval)
{
    OGL_ViewerRotate(viewer, newval - oldval, ROTATE_Y, FALSE, TRUE);
}

static void OGL_ViewerZScrollProc(Nlm_BaR sb, Nlm_GraphiC group, Nlm_Int2 newval, Nlm_Int2 oldval)
{
    Nlm_SlatE viewer = (Nlm_SlatE)Nlm_GetObjectExtra( sb );
    OGL_ViewerRotate(viewer, newval - oldval, ROTATE_Z, FALSE, TRUE);
}




/*
 *  MOUSE EVENT HANDLERS
 */



/*
 * MOVE
 */

static void OGL_Move_DrawTrace(MA_TracePtr trace)
{
    if ( Nlm_EqualPt(trace->start, trace->end) )
        return;
    
#ifdef WIN_MOTIF
    Nlm_SetColor( 0xf1 );
#endif
    Nlm_InvertMode();  /* turn on xor */
    Nlm_DrawLine(trace->start, trace->end);  /* draw to current hdc */
    Nlm_CopyMode();  /* turn off xor */
}


static void OGL_Move_PressMA(MAPtr ma,
                             MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    trace->start = trace->end = point;
}


static void OGL_Move_DragMA(MAPtr ma,
                            MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Move_DrawTrace( trace );
    trace->end = point;
    OGL_Move_DrawTrace( trace );
}


static void OGL_Move_ReleaseMA(MAPtr ma,
                               MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Move_DrawTrace( trace );
    trace->end = point;
    
    if ( Nlm_EqualPt(trace->start, trace->end) )
        return;
    
    { /* do the move transform */
        TOGL_Data * OGL_Data = MAToOGL( ma );
        
        OGL_Move3D(OGL_Data, (Int2)(trace->end.x - trace->start.x),
            (Int2)(trace->end.y - trace->start.y));
        OGL_DrawViewer3D(OGL_Data);
        
    }
    
    trace->start = trace->end;
}


static void OLG_Move_CancelMA(MAPtr ma,
                              MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Move_DrawTrace( trace );
}


/*
 * ZOOM
 */

static void OGL_Zoom_DrawTrace(MA_TracePtr trace)
{
    Nlm_RecT rubber_box;
    if ( Nlm_EqualPt(trace->start, trace->end) )
        return;
    
#ifdef WIN_MOTIF
    Nlm_SetColor( 0xf1 );
#endif
    Nlm_InvertMode();
    Nlm_LoadRect(&rubber_box, trace->start.x, trace->start.y,
        trace->end.x,   trace->end.y);
    Nlm_FrameRect( &rubber_box );  /* draw the frame */
    Nlm_CopyMode();
}

static void OGL_Zoom_PressMA(MAPtr ma,
                             MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    trace->start = trace->end = point;
}


static void OGL_Zoom_DragMA(MAPtr ma,
                            MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Zoom_DrawTrace( trace );
    trace->end = point;
    OGL_Zoom_DrawTrace( trace );
}


static void OGL_Zoom_ReleaseMA(MAPtr ma,
                               MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Zoom_DrawTrace( trace );
    trace->end = point;
    
    if ( Nlm_EqualPt(trace->start, trace->end) )
        return;
    
    { /* do the zoom */
        TOGL_Data * OGL_Data = MAToOGL( ma );
        
        OGL_Zoom3D(OGL_Data, trace->start.x, trace->start.y,
            trace->end.x,   trace->end.y);
        OGL_DrawViewer3D(OGL_Data);
        
    }
}

static void OGL_Zoom_CancelMA(MAPtr ma,
                              MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    OGL_Zoom_DrawTrace( trace );
}


/*
 * ROTATE
 */


static void OGL_Rotate_PressMA(MAPtr ma,
                               MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    trace->start = point;
}



static void OGL_Rotate_DragMA(MAPtr ma,
                              MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    TOGL_Data * OGL_Data;

    OGL_RotatePivots3DPtr pivot;
    if ( Nlm_EqualPt(trace->start, point) )
        return;
    
    OGL_Data = MAToOGL( ma );
    pivot = (OGL_RotatePivots3DPtr)extra;
    
    OGL_ViewerRotate((Nlm_SlatE)OGL_Data->Panel,
        (Int4)((180.0 * (point.x - trace->start.x)) / OGL_Data->Width),
        pivot->H, TRUE, FALSE);
    
    OGL_ViewerRotate((Nlm_SlatE)OGL_Data->Panel,
        (Int4)((180.0 * (trace->start.y - point.y)) / OGL_Data->Height),
        pivot->V, TRUE, TRUE);
    trace->start = point;
}


/*
 * RESET
 */

static void OGL_ResetMA(MAPtr ma,
                        MA_TracePtr trace, Nlm_PoinT point, Nlm_VoidPtr extra)
{
    VERIFY ( MA_UnsetAll( ma ) );
}


static Nlm_Boolean OGL_SetStdMouse(TOGL_Data * OGL_Data, Nlm_enumStdMAOGL action)
{
    if(OGL_Data == NULL) return FALSE;
    return MA_SetGroup( OGL_Data->ma_std_group[action] );
}



/* Initialize MA for the viewer
*/

static Nlm_Boolean OGL_InitializeMA(TOGL_Data * OGL_Data)
{
    MAPtr ma = OGL_Data->ma; 
    
    /* rotate */
    MActionPtr rotate_press   =
        MA_AddAction(ma, MK_Normal, MA_Press,   OGL_Rotate_PressMA, NULL, NULL);
    
    static OGL_RotatePivots3D RotateDrag_YX = { ROTATE_Y, ROTATE_X };
    MActionPtr rotate_drag_YX =
        MA_AddAction(ma, MK_Normal, MA_Drag, OGL_Rotate_DragMA, &RotateDrag_YX, NULL);
    MA_GroupPtr rotate_group_YX = MA_AddGroup(ma, "Rotate_YX",
        rotate_press,   MA_ONLY,
        rotate_drag_YX, MA_ONLY,
        NULL);
    
    static OGL_RotatePivots3D RotateDrag_ZX = { ROTATE_Z, ROTATE_X };
    MActionPtr rotate_drag_ZX =
        MA_AddAction(ma, MK_Normal, MA_Drag, OGL_Rotate_DragMA, &RotateDrag_ZX, NULL);
    MA_GroupPtr rotate_group_ZX = MA_AddGroup(ma, "Rotate_ZX",
        rotate_press,   MA_ONLY,
        rotate_drag_ZX, MA_ONLY,
        NULL);
    
    static OGL_RotatePivots3D RotateDrag_YZ = { ROTATE_Y, ROTATE_Z };
    MActionPtr rotate_drag_YZ =
        MA_AddAction(ma, MK_Normal, MA_Drag, OGL_Rotate_DragMA, &RotateDrag_YZ, NULL);
    MA_GroupPtr rotate_group_YZ = MA_AddGroup(ma, "Rotate_YZ",
        rotate_press,   MA_ONLY,
        rotate_drag_YZ, MA_ONLY,
        NULL);
    
    /* move */
    MActionPtr move_press   =
        MA_AddAction(ma, MK_Shift, MA_Press,   OGL_Move_PressMA,   NULL, NULL);
    MActionPtr move_drag    =
        MA_AddAction(ma, MK_Shift, MA_Drag,    OGL_Move_DragMA,    NULL, NULL);
    MActionPtr move_release =
        MA_AddAction(ma, MK_Shift, MA_Release, OGL_Move_ReleaseMA, NULL, NULL);
    MActionPtr move_cancel  =
        MA_AddAction(ma, MK_Shift, MA_Cancel,  OLG_Move_CancelMA,  NULL, NULL);
    
    MA_GroupPtr move_group = MA_AddGroup(ma, "Move",
        move_press,   MA_ONLY,
        move_drag,    MA_ONLY,
        move_release, MA_ONLY,
        move_cancel,  MA_ONLY,
        NULL);
    
    /* zoom */
    MActionPtr zoom_press   =
        MA_AddAction(ma, MK_Ctrl, MA_Press,   OGL_Zoom_PressMA,   NULL, NULL);
    MActionPtr zoom_drag    =
        MA_AddAction(ma, MK_Ctrl, MA_Drag,    OGL_Zoom_DragMA,    NULL, NULL);
    MActionPtr zoom_release =
        MA_AddAction(ma, MK_Ctrl, MA_Release, OGL_Zoom_ReleaseMA, NULL, NULL);
    MActionPtr zoom_cancel  =
        MA_AddAction(ma, MK_Ctrl, MA_Cancel,  OGL_Zoom_CancelMA,  NULL, NULL);
    
    MA_GroupPtr zoom_group = MA_AddGroup(ma, "Zoom",
        zoom_press,   MA_ONLY,
        zoom_drag,    MA_ONLY,
        zoom_release, MA_ONLY,
        zoom_cancel,  MA_ONLY,
        NULL);
    
    /* miscellaneous actions */
/*
 this is done in the main program.  move it here after deleting viewer3d
    MActionPtr bg_hl_dclick =
        MA_AddAction(ma, MK_Normal, MA_DClick,  NULL, NULL,
        "Highlight-Prim or Background");
*/
    
    /* this group disables all mouse actions when set */
    MActionPtr reset_init =
        MA_AddAction(ma, MK_Normal, MA_Init,  OGL_ResetMA, NULL, NULL);
    
    MA_GroupPtr reset_group = MA_AddGroup(ma, "No Action",
        reset_init, MA_SHARED,
        NULL);
    
    if(OGL_Data == NULL) return FALSE;
    
    {{ /* "No-Action"s */
        int i, j;
        for (i = 0;  i < MK_Default;  i++)
            for (j = 0;  j < MA_Init;     j++)
            {
                VERIFY ( MA_AddAction(ma, (enumMKey)i, (enumMAction)j,
                    DoNothingMA,   NULL, "No Action") );
            }
    }}
    
    /* register the set of standard 3D-viewer groups */
    OGL_Data->ma_std_group[MouseOGL_DoNothing] = reset_group;
    OGL_Data->ma_std_group[MouseOGL_RotateYX ] = rotate_group_YX;
    OGL_Data->ma_std_group[MouseOGL_RotateZX ] = rotate_group_ZX;
    OGL_Data->ma_std_group[MouseOGL_RotateYZ ] = rotate_group_YZ;
    OGL_Data->ma_std_group[MouseOGL_Move     ] = move_group;
    OGL_Data->ma_std_group[MouseOGL_Zoom     ] = zoom_group;
    
    /* Test, Setup defaults and Link viewer panel to MA */
    if (!rotate_press  ||
        !rotate_drag_YX  ||  !rotate_group_YX  ||
        !rotate_drag_ZX  ||  !rotate_group_ZX  ||
        !rotate_drag_YZ  ||  !rotate_group_YZ  ||
        !move_press      ||  !move_drag        ||  !move_release  ||
        !move_cancel     ||  !move_group       ||
        !zoom_press      ||  !zoom_drag        ||  !zoom_release  ||
        !zoom_cancel     ||  !zoom_group       ||
/*        !bg_hl_dclick    ||*/
        !reset_group     ||
        
        !OGL_SetStdMouse(OGL_Data, MouseOGL_RotateYX)  ||
        !OGL_SetStdMouse(OGL_Data, MouseOGL_Move    )  ||
        !OGL_SetStdMouse(OGL_Data, MouseOGL_Zoom    )  ||
/*        !MA_SetAction(bg_hl_dclick, FALSE)  ||*/
        !MA_LinkPanel(ma, OGL_Data->Panel))
    {
        MA_Reset( ma );
        return FALSE;
    }
    
    return TRUE;
}


/*
 *  Doing selection in OpenGL
 */

void OGL_Select(TOGL_Data * OGL_Data, Nlm_Boolean SelectMode)
{
    if (OGL_Data == NULL) return;
    OGL_Data->SelectMode = SelectMode;
    return;
}


void OGL_LoadName(Nlm_VoidPtr PtrValue)
/* load a pointer onto the name stack.  compensate for possible long long */
{
    Nlm_Int4 i;
    
    for(i = 0; i < sizeof(Nlm_VoidPtr)/sizeof(GLuint); i++ ) glPopName();

    for(i = 0; i < sizeof(Nlm_VoidPtr)/sizeof(GLuint); i++ ) 
        glPushName((GLuint)(((long)PtrValue) >> (i * sizeof(GLuint) * 8)));  /* 64 bits? */
}


Nlm_VoidPtr OGL_Hit(TOGL_Data * OGL_Data)
/* this function looks through the hit stack and extracts the nearest hit */
{
    GLuint * Hits, Names, i = 0, j, ZMax = 0;
    long ZMaxName = 0;  /* this is a hack, but should work for 64 bits */

    if(OGL_Data == NULL || OGL_Data->SelectBuffer == NULL) return NULL;

    Hits = (GLuint *) OGL_Data->SelectBuffer;
    while (i < OGL_Data->SelectHits) {
        Names = Hits[i];
        i += 2; /* skip to z max */
        if( Hits[i] > ZMax) {
            ZMax = Hits[i++];
            for(j = 0; j < sizeof(Nlm_VoidPtr)/sizeof(GLuint); j++ ) {
                ZMaxName << (sizeof(GLuint) * 8);
                ZMaxName |= Hits[i+j];
            }
            i += Names * sizeof(Nlm_VoidPtr)/sizeof(GLuint);  /* currently only looks at the top of the name stack */
        }
        else i += Names * sizeof(Nlm_VoidPtr)/sizeof(GLuint) + 1;
    }
return (Nlm_VoidPtr) ZMaxName;
}
       

void OGL_SetSelectPoint(TOGL_Data * OGL_Data, Nlm_PoinT Point)
{
    if(OGL_Data == NULL) return;
    OGL_Data->SelectPoint.x = Point.x;
    OGL_Data->SelectPoint.y = Point.y;
}


/*
 *  functions used to create, do, and finish drawing
 */


static void OGL_DeleteViewer3D(TOGL_Data * OGL_Data)
/* delete OGL_Data */
{
    if (OGL_Data == NULL)
        return;
    
    MA_Destroy( OGL_Data->ma );
    
    /* to do: someone else is getting rid of Panel
    if ( OGL_Data->Panel )
    Nlm_Remove( OGL_Data->Panel );
    */

    /* delete the display lists */
    glDeleteLists(OGL_Data->Layers->FirstLayer - 1, OGLMAXLAYERS);
    
    /* free items on the heap */
    MemFree( OGL_Data->Layers );
    MemFree( OGL_Data->ModelMatrix );
    MemFree( OGL_Data );
}


static void OGL_ResetViewerProc_CB(Nlm_PaneL panel)
{
    TOGL_Data * OGL_Data = PanelToOGL( panel );
    if ( !OGL_Data  )
        return;
    
    OGL_Data = NULL;
    OGL_DeleteViewer3D( PanelToOGL( panel ) );
}


TOGL_Data * OGL_CreateViewer(Nlm_GrouP prnt,
                             Uint2Ptr width, Uint2 height,
                             Int4 flags,
                             Nlm_MenU ma_group_menu, Nlm_MenU ma_action_menu,
                             Nlm_MAInitOGLFunc ma_init_func,
                             VoidPtr          ma_init_data) 
                             /* initialize the OpenGL library */
{
    TOGL_Data * OGL_Data;
    Nlm_Uint2 x_width;    
    
    OGL_Data = (TOGL_Data *)MemNew( sizeof(TOGL_Data) );
    if (OGL_Data == NULL) return NULL;
    
    OGL_Data->ModelMatrix = (Nlm_VoidPtr) MemNew( 16 * sizeof(GLdouble));
    if (OGL_Data->ModelMatrix == NULL) return NULL;
    
    OGL_Data->Layers = (TOGL_Layers *) MemNew( sizeof(TOGL_Layers));
    if (OGL_Data->Layers == NULL) return NULL;
    
    OGL_Data->Layers->FirstLayer = glGenLists((GLsizei)OGLMAXLAYERS);
    OGL_Data->Layers->FirstLayer++;  /* avoid weird bug in windows OpenGL */
    OGL_Data->Layers->LastLayer = OGL_Data->Layers->FirstLayer;
    OGL_Data->Layers->SelectedLayer = 0;
    OGL_SetLayers(OGL_Data, FALSE);  /* null all the layers out */
    
    OGL_Data->IsPlaying = FALSE;  /* animation off */
    OGL_Data->ParentWindow = Nlm_ParentWindow( (Nlm_Handle)prnt );
    OGL_Data->PaletteExpanded = NULL;
    OGL_Data->PaletteIndex = NULL;
    OGL_Data->NewPalette = FALSE;

    OGL_Data->SelectMode = FALSE;
    OGL_Data->SelectBuffer = (Nlm_VoidPtr) MemNew(OGLSELECTBUFFER * sizeof(GLuint));
    if (OGL_Data->SelectBuffer == NULL) return NULL;
    
    OGL_Data->Panel = Nlm_Autonomous3DPanel(prnt,
        (Int2)*width, (Int2)height,
        OGL_DrawViewer3D_CB,
        ((flags & Y_ROTATE_SBAR) ?
OGL_ViewerVScrollProc : NULL),
                        ((flags & X_ROTATE_SBAR) ?
OGL_ViewerHScrollProc : NULL),
                        sizeof(MAPtr), OGL_ResetViewerProc_CB, NULL, &OGL_Data->IndexMode);
 
    
    if (flags & Z_ROTATE_SBAR)
    {
        OGL_Data->Z_rotate = Nlm_ScrollBar(prnt, 1, 0, OGL_ViewerZScrollProc);
        if ( !OGL_Data->Z_rotate ) {
            MemFree(OGL_Data);
            if ( OGL_Data->Panel )
                Nlm_Remove( OGL_Data->Panel );
            return NULL;
        }
        Nlm_SetObjectExtra(OGL_Data->Z_rotate, OGL_Data->Panel, NULL);
    }
    
    {
        Nlm_RecT rect;
        Nlm_GetPosition(OGL_Data->Panel, &rect);
        rect.right  = (Int2)(rect.left + *width);
        rect.bottom = (Int2)(rect.top  + height);
        OGL_SetPosition3D(OGL_Data, &rect);
        x_width = (Uint2)(rect.right - rect.left);
    }
    
    if (flags & X_ROTATE_SBAR) {
        Nlm_BaR sb = Nlm_GetSlateHScrollBar( (Nlm_SlatE)OGL_Data->Panel );
        Nlm_CorrectBarValue(sb, 0);
        Nlm_SetRange(sb, 10, 10, 360);
    }
    if (flags & Y_ROTATE_SBAR) {
        Nlm_BaR sb = Nlm_GetSlateVScrollBar( (Nlm_SlatE)OGL_Data->Panel );
        Nlm_CorrectBarValue(sb, 0);
        Nlm_SetRange(sb, 10, 10, 360);
    }
    if (flags & Z_ROTATE_SBAR) {
        Nlm_SetRange(OGL_Data->Z_rotate, 10, 10, 360);
        Nlm_CorrectBarValue(OGL_Data->Z_rotate, 180);
    }
    
    OGL_Data->ma = MA_Create(ma_group_menu, ma_action_menu);
    MA_SetExtra(OGL_Data->ma, OGL_Data);
    
    if ( !OGL_InitializeMA( OGL_Data ) ) {
        MemFree(OGL_Data);
        if ( OGL_Data->Z_rotate ) Nlm_Remove( OGL_Data->Z_rotate );
        if ( OGL_Data->Panel ) Nlm_Remove( OGL_Data->Panel );
        return NULL;
    }
    
    if (ma_init_func  &&  !(*ma_init_func)(OGL_Data->ma, ma_init_data)) {
        MemFree(OGL_Data);
        if ( OGL_Data->Z_rotate ) Nlm_Remove( OGL_Data->Z_rotate );
        if ( OGL_Data->Panel ) Nlm_Remove( OGL_Data->Panel );
        return NULL;
    }
    
    *width = x_width;
    
    
    
    /* set up the font */
#ifdef WIN32
    {
        HDC hdc; 
        SIZE TextSize;
        
        hdc = wglGetCurrentDC();
        SelectObject (hdc, GetStockObject (SYSTEM_FONT));  
        wglUseFontBitmaps (hdc, 0, 255, OGLFONTBASE);
        /* todo: explore using outlines if cross platform.  The problem is that it isn't supported on X */
        GetTextExtentPoint32(hdc, "A", 1, &TextSize);  /* Get the size of an average character */
        
        OGL_Data->SpaceWidth = TextSize.cx;
        OGL_Data->SpaceHeight = TextSize.cy;
        OGL_Data->Space = (Nlm_VoidPtr) MemNew(TextSize.cy * ((TextSize.cx - 1)%8 + 1));
        MemSet(OGL_Data->Space, 0, (size_t) (TextSize.cy * ((TextSize.cx - 1)%8 + 1)));
    }
#endif
    
    return OGL_Data;
}


void OGL_DrawViewer3D(TOGL_Data * OGL_Data)
/* does the drawing */
{
    GLfloat LightPosition[4];
    GLfloat Ambient[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat Diffuse[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat Specular[4] = { 0.3f, 0.3f, 0.3f, 0.3f };
    GLfloat Background[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    Nlm_Uint1 * red, * green, * blue;
    Nlm_Int4 TotalColors, i;
    Nlm_Uint4 iList;
    ValNodePtr Palette, PaletteIndex; 
    GLint Viewport[4];
    
    if(OGL_Data == NULL) return;
    if(OGL_Data->SelectMode) {
        glSelectBuffer((GLsizei)OGLSELECTBUFFER, (GLuint *) OGL_Data->SelectBuffer);
        glRenderMode(GL_SELECT);
        glGetIntegerv(GL_VIEWPORT, Viewport);
        glInitNames();
        for(i = 0; i < sizeof(Nlm_VoidPtr)/sizeof(GLuint); i++ ) glPushName(0);
    }

    
    /* set up the lighting */

    if(OGL_Data->IndexMode) {
        TotalColors = ValNodeLen(OGL_Data->PaletteExpanded);
        if(TotalColors && OGL_Data->NewPalette) {
            OGL_Data->NewPalette = FALSE;
            red = Calloc((size_t)TotalColors + 1, sizeof(Uint1));
            blue = Calloc((size_t)TotalColors + 1, sizeof(Uint1));
            green = Calloc((size_t)TotalColors + 1, sizeof(Uint1));
            
            red[0] = OGL_Data->Background.rgb[0];
            green[0] = OGL_Data->Background.rgb[1];
            blue[0] = OGL_Data->Background.rgb[2];
            
            for(Palette = OGL_Data->PaletteExpanded, i = 1; Palette; Palette = Palette->next, i++) {
                red[i] = ((TOGL_ColorCell *)(Palette->data.ptrvalue))->rgb[0];
                green[i] = ((TOGL_ColorCell *)(Palette->data.ptrvalue))->rgb[1];
                blue[i] = ((TOGL_ColorCell *)(Palette->data.ptrvalue))->rgb[2];
            }
            
            Nlm_Set3DColorMap (OGL_Data->Panel, (Uint1)(TotalColors + 1), red, green, blue);

            /* increment the palette index to allow for the background */
            for(PaletteIndex = OGL_Data->PaletteIndex; PaletteIndex; PaletteIndex = PaletteIndex->next) {
                ((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->Begin += 1;
                ((TOGL_PaletteIndex *)(PaletteIndex->data.ptrvalue))->End += 1;
            }
           
            MemFree(red);
            MemFree(green);
            MemFree(blue);
        }
        glClearIndex(0.0);
    }
    else {
        ColorCell2Array(Background, &OGL_Data->Background);
        glClearColor( Background[0], Background[1], Background[2], Background[3]);
    }
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);  /* turn on one light, at infinity */
    glEnable(GL_DEPTH_TEST);
    /* turn on ambient light for lines */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, Diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, Specular);
    
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (OGL_Data->SelectMode) gluPickMatrix((GLdouble)(OGL_Data->SelectPoint.x - 5),
        (GLdouble)(Viewport[3] - OGL_Data->SelectPoint.y + 25), 1.0, 1.0, Viewport);
    {
        GLfloat Height = (GLfloat)(OGL_Data->ViewPoint[3]/18.0), Width =
            (GLfloat)(OGL_Data->ViewPoint[3]/18.0);
        
        if(OGL_Data->Width > OGL_Data->Height) Width *= (GLfloat)OGL_Data->Width/OGL_Data->Height;
        else Height *= (GLfloat)OGL_Data->Height/OGL_Data->Width;
        
        
        glFrustum(-Width, Width, -Height, Height, OGL_Data->ViewPoint[3]/3.0, 5*OGL_Data->MaxSize);  
    }
    glMatrixMode(GL_MODELVIEW);
    
    /* set up the lighting */
    LightPosition[0] = 0.0;
    LightPosition[1] = 0.0;
    LightPosition[2] = (GLfloat)3.0*OGL_Data->MaxSize;
    LightPosition[3] = 1.0;
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
    gluLookAt(OGL_Data->ViewPoint[0], OGL_Data->ViewPoint[1], 3.0*OGL_Data->MaxSize,
        OGL_Data->ViewPoint[0], OGL_Data->ViewPoint[1], 0.0,0.0,1.0,0.0); 
    glMultMatrixd((GLdouble *)OGL_Data->ModelMatrix);
    

    /* selectively display the layers */
    if(OGL_Data->Layers->SelectedLayer) {
        glCallList(OGL_Data->Layers->SelectedLayer);
    }
    else {
        for(iList = OGL_Data->Layers->FirstLayer;
        iList <= OGL_Data->Layers->LastLayer; iList++) {
            if(OGL_GetLayer(OGL_Data, iList)) glCallList(iList);
        }
    }
    
    glFlush();

    if(OGL_Data->SelectMode) OGL_Data->SelectHits = glRenderMode(GL_RENDER);
    else {
    /* swap when double buffering */
#ifdef WIN32
        SwapBuffers(wglGetCurrentDC());
#endif
    }
    
}


/*
*  Helper functions for drawing
*/

void OGL_Redraw(TOGL_Data * OGL_Data)
{
      Nlm_RecT   r;
      Nlm_WindoW tempPort;

      if(OGL_Data == NULL) return;
      tempPort = Nlm_SavePort( OGL_Data->Panel );
      Nlm_Select( OGL_Data->Panel );
      Nlm_ObjectRect(OGL_Data->Panel, &r);
      Nlm_InvalRect( &r );
      Nlm_RestorePort( tempPort );
      return;
}


void OGL_DataFromBound(TOGL_Data * OGL_Data)
/* create origin and bound size from the bound box */
{
    if(OGL_Data == NULL) return;
    OGL_Data->Translate[0] = (Nlm_FloatLo)((OGL_Data->BoundBox.x[0] + OGL_Data->BoundBox.x[1])/2.0);
    OGL_Data->Translate[1] = (Nlm_FloatLo)((OGL_Data->BoundBox.y[0] + OGL_Data->BoundBox.y[1])/2.0);
    OGL_Data->Translate[2] = (Nlm_FloatLo)((OGL_Data->BoundBox.z[0] + OGL_Data->BoundBox.z[1])/2.0);
    OGL_Data->MaxSize = (Nlm_FloatLo)fabs(OGL_Data->BoundBox.x[0] - OGL_Data->BoundBox.x[1]);
    if( fabs(OGL_Data->BoundBox.y[0] - OGL_Data->BoundBox.y[1]) > OGL_Data->MaxSize) 
        OGL_Data->MaxSize = (FloatLo)fabs(OGL_Data->BoundBox.y[0] - OGL_Data->BoundBox.y[1]);
    if( fabs(OGL_Data->BoundBox.z[0] - OGL_Data->BoundBox.z[1]) > OGL_Data->MaxSize) 
        OGL_Data->MaxSize = (FloatLo)fabs(OGL_Data->BoundBox.z[0] - OGL_Data->BoundBox.z[1]);
}


void OGL_Reset(TOGL_Data * OGL_Data)
/* reset the transform matrix */
{
    if(OGL_Data == NULL) return;
    OGL_DataFromBound(OGL_Data);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    /* set up the initial view point */
    OGL_Data->ViewPoint[0] = 0.0;
    OGL_Data->ViewPoint[1] = 0.0;
    OGL_Data->ViewPoint[2] = (GLfloat)3.0*OGL_Data->MaxSize;
    
    
    glTranslatef(-(GLfloat)(OGL_Data->Translate[0]), -(GLfloat)(OGL_Data->Translate[1]), -(GLfloat)(OGL_Data->Translate[2]));
    
    glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *)OGL_Data->ModelMatrix); 
}


Boolean OGL_SetPosition3D(TOGL_Data * OGL_Data, Nlm_RectPtr rect)
/* resizes and positions the 3D window */
{
    Nlm_Int4 width  = rect->right - rect->left;
    Nlm_Int4 height = rect->bottom - rect->top;

    if(OGL_Data == NULL || rect == NULL) return FALSE;
    
    if ( OGL_Data->Z_rotate )
    {
        rect->top += Nlm_hScrollBarHeight;
        height    -= Nlm_hScrollBarHeight;
    }
    
    if ( Nlm_GetSlateVScrollBar((Nlm_SlatE)OGL_Data->Panel) )
        width  -= Nlm_vScrollBarWidth  + 3;
    if ( Nlm_GetSlateHScrollBar((Nlm_SlatE)OGL_Data->Panel) )
        height -= Nlm_hScrollBarHeight + 3;
    
    if (width < 16  ||  height < 16 ) return FALSE;
    
    OGL_Data->Width  = (Nlm_Int2)width;
    OGL_Data->Height = (Nlm_Int2)height;
    
    glViewport(0, 0, OGL_Data->Width + 4, OGL_Data->Height + 4);
    
    rect->right  = (Nlm_Int2)(rect->left + OGL_Data->Width  + 3);
    rect->bottom = (Nlm_Int2)(rect->top  + OGL_Data->Height + 3);
    if ( Nlm_GetSlateVScrollBar((Nlm_SlatE)OGL_Data->Panel) )
        rect->right  += Nlm_vScrollBarWidth;
    if ( Nlm_GetSlateHScrollBar((Nlm_SlatE)OGL_Data->Panel) )
        rect->bottom += Nlm_hScrollBarHeight;
    
    Nlm_SetPosition(OGL_Data->Panel, rect);  /* resize the panel */
    Nlm_ProcessUpdatesFirst( FALSE );
    Nlm_AdjustPrnt(OGL_Data->Panel, rect, FALSE);
    
    if ( OGL_Data->Z_rotate )  /* if there is a z rotation scroll bar */
    {
        rect->bottom = rect->top;
        rect->top   -= Nlm_hScrollBarHeight;
        Nlm_SetPosition(OGL_Data->Z_rotate, rect);
    }
    
    return TRUE;
}


void OGL_ClearOGL_Data (TOGL_Data * OGL_Data)
/* clear the transforms and bound box in OGL_Data structure */
{
    Nlm_Int4 i;

    if(OGL_Data == NULL) return;
    for (i=0; i<3; i++) {
        OGL_Data->Translate[i] = 0.0;
    }
    OGL_ClearBoundBox(&(OGL_Data->BoundBox));
    OGL_Data->MaxSize = 2.0*OGL_DEFAULT_SIZE;
    OGL_Data->Background.rgb[0] = 0;
    OGL_Data->Background.rgb[1] = 0;
    OGL_Data->Background.rgb[2] = 0;
}


void OGL_ClearBoundBox(TOGL_BoundBox * BoundBox)
/* initialize a bounds box */
{
    Nlm_Int4 i;
    
    if (BoundBox == NULL) return;
    for (i=0; i<2; i++ ) {
        BoundBox->x[i] = (Nlm_FloatLo)((i*2-1)*OGL_DEFAULT_SIZE);
        BoundBox->y[i] = (Nlm_FloatLo)((i*2-1)*OGL_DEFAULT_SIZE);
        BoundBox->z[i] = (Nlm_FloatLo)((i*2-1)*OGL_DEFAULT_SIZE);
    }        
    BoundBox->set = TRUE;  /* we've set up the default values */
}



#endif /* _OPENGL */
