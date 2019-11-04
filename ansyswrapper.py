import os
import datetime
import tempfile

from subprocess import call
import numpy as np


class ansyswrapper:
    __matid = 0
    __csid = 10

    def __init__(self, ans_hi_version=250, ans_low_version=50, anslic='ane3fl', infile='input.dat',
                 outfile='output.dat', projdir=tempfile.gettempdir(), isBatch=True, jobname=None):

        self.ans_hi_version = ans_hi_version
        self.ans_low_version = ans_low_version

        self.anslic = anslic
        self.inputfile = infile
        self.outputfile = outfile

        if not jobname:
            dt = datetime.datetime.now()
            self.jobname = 'jobname{0}-{1}-{2}--{3}-{4}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute)
        else:
            self.jobname = jobname

        self.apdl = ""
        self.apdl += "FINISH\n"
        self.apdl += "/CLEAR,START\n"
        self.apdl += "/prep7\n"
        self.apdl += "it_num = 1\n"
        self.apdl += self.__post_new
        self.projdir = projdir
        self.__isBatch = isBatch

        self.apdl += """
            !--------------- randomize ----------------
            *GET,DIM,ACTIVE,0,TIME,WALL
            DIM=DIM*3600
            *DIM,DUMMY,ARRAY,DIM
            *VFILL,DUMMY(1),RAND
            *DEL,DIM
            *DEL,DUMMY
            !--------------- randomize ----------------\n"""

        self.apdl += "*DMAT, mat_s, D, Alloc, 3, 3, incore\n"
        self.apdl += "*VEC, vec_b, D, alloc, 3\n"
        self.apdl += "*VEC, vec_x, D, alloc, 3\n"

    def saveToFile(self, filename):
        f = open(filename, mode='w')
        f.write(self.apdl)
        f.close()

    def getNP(self):
        return 4  # os.environ['NUMBER_OF_PROCESSORS']

    def findPathVersion(self):

        __ansversion = -1
        path = ""
        for i in range(self.ans_hi_version, self.ans_low_version, -1):
            if os.environ.get('ANSYS{0}_DIR'.format(i)):
                path = os.environ.get('ANSYS{0}_DIR'.format(i))
                __ansversion = i
                break

        if __ansversion == -1:
            print("Ansys not found")
            return None

        path += '\\bin\\' + os.environ['ANSYS_SYSDIR'] + '\\ANSYS{0}'.format(__ansversion)

        return path

    def defaultArgs(self):
        if self.__isBatch:
            return '-b -p {0} -np {1} -dir {2} -j {3} -s noread -i {4} -o {5} -d win32'.format(
                self.anslic, self.getNP(), self.projdir, self.jobname, self.inputfile, self.outputfile
            )
        else:
            return '-g -p {0} -np {1} -dir {2} -j {2} -s read -d win32'.format(
                self.anslic, self.getNP(), self.projdir, self.jobname
            )

    def run(self, apdl=None):

        self.apdl += "/EXIT, NOSAVE,\n"
        cwd = os.getcwd()
        os.chdir(self.projdir)
        self.saveToFile(self.projdir + "\\" + self.inputfile)

        retcode = call([self.findPathVersion()] + self.defaultArgs().split(" "))
        os.chdir(cwd)

        exitcodes = dict()

        # exitcodes[0] = 'Normal Exit'
        exitcodes[1] = 'Stack Error'
        exitcodes[2] = 'Stack Error'
        exitcodes[3] = 'Stack Error'
        exitcodes[4] = 'Stack Error'
        exitcodes[5] = 'Command Line Argument Error'
        exitcodes[6] = 'Accounting File Error'
        exitcodes[7] = 'Auth File Verification Error'
        exitcodes[8] = 'Error in ANSYS or End-of-run'
        exitcodes[11] = 'User Routine Error'
        exitcodes[12] = 'Macro STOP Command'
        exitcodes[14] = 'XOX Error'
        exitcodes[15] = 'Fatal Error'
        exitcodes[16] = 'Possible Full Disk'
        exitcodes[17] = 'Possible Corrupted or Missing File'
        exitcodes[18] = 'Possible Corrupted DB File'
        exitcodes[21] = 'Authorized Code Section Entered'
        exitcodes[25] = 'Unable to Open X11 Server'
        exitcodes[30] = 'Quit Signal'
        exitcodes[31] = 'Failure to Get Signal'
        exitcodes[32] = 'System-dependent Error'

        if retcode > 32:
            exitcodes[retcode] = 'Unknown Error. Check for *.lock files in working directory and delete it'

        if retcode in exitcodes:
            print('------ANSYS ERROR EXIT CODE-------')
            print('Ansys exit code = {0}, with message: {1}'.format(retcode, exitcodes[retcode]))
            print('----------------------------------')
            print('Terminating.......')
            exit(retcode)
        return retcode

    def rectangle(self, x1, y1, x2, y2):
        self.apdl += "RECTNG,{0},{1},{2},{3},\n".format(x1, x2, y1, y2)

    def circle(self, x, y, rad):
        self.apdl += "CYL4, {0}, {1}, {2}\n".format(x, y, rad)

    def ellipse(self, x, y, r1, r2):
        self.apdl += "ASEL, NONE\n"
        self.circle(x, y, 1)
        self.apdl += "ARSCALE, ALL, , , {0}, {1}, , , 1, 1\n".format(r1, r2)
        self.apdl += "ASEL, ALL\n"

    def setFEByNum(self, num):
        self.apdl += "ET, 1, {0}\n".format(num)

        # self.apdl += "KEYOPT, 1, 3, 2\n"

    def createIsotropicMat(self, E, nu):
        self.__matid += 1
        self.apdl += "MPTEMP,, , , , , , ,\n"
        self.apdl += "MPTEMP, 1, 0\n"
        self.apdl += "MPDATA, EX, {1},, {0}\n".format(E, self.__matid)
        self.apdl += "MPDATA, PRXY, {1},, {0}\n".format(nu, self.__matid)
        return self.__matid

    def overlapAreas(self):
        self.apdl += "AOVLAP,ALL\n"

    def createOrtotropicMat(self, c11, c12, c13, c22, c23, c33, c44, c55=None, c66=None):
        self.__matid += 1

        if c55 == None:
            c55 = c44

        if c66 == None:
            c66 = c44

        Ex = (c11 * c22 * c33 + 2 * c23 * c12 * c13 - c11 * c23 ** 2 - c22 * c13 ** 2 - c33 * c12 ** 2) / (
                c22 * c33 - c23 ** 2)
        Ey = (c11 * c22 * c33 + 2 * c23 * c12 * c13 - c11 * c23 ** 2 - c22 * c13 ** 2 - c33 * c12 ** 2) / (
                c11 * c33 - c13 ** 2)
        Ez = (c11 * c22 * c33 + 2 * c23 * c12 * c13 - c11 * c23 ** 2 - c22 * c13 ** 2 - c33 * c12 ** 2) / (
                c11 * c22 - c12 ** 2)
        nuxy = (c12 * c33 - c13 * c23) / (c22 * c33 - c23 ** 2)
        nuxz = (c22 * c13 - c12 * c23) / (c22 * c33 - c23 ** 2)
        nuyz = (c11 * c23 - c12 * c13) / (c11 * c33 - c13 ** 2)

        self.apdl += "MPTEMP,, , , , , , ,\n"
        self.apdl += "MPTEMP, 1, 0\n"
        self.apdl += "MPDATA, EX, {1},, {0}\n".format(Ex, self.__matid)
        self.apdl += "MPDATA, EY, {1},, {0}\n".format(Ey, self.__matid)
        self.apdl += "MPDATA, EZ, {1},, {0}\n".format(Ez, self.__matid)
        self.apdl += "MPDATA, PRXY, {1},, {0}\n".format(nuxy, self.__matid)
        self.apdl += "MPDATA, PRYZ, {1},, {0}\n".format(nuyz, self.__matid)
        self.apdl += "MPDATA, PRXZ, {1},, {0}\n".format(nuxz, self.__matid)
        self.apdl += "MPDATA, GXY, {1},, {0}\n".format(c44, self.__matid)
        self.apdl += "MPDATA, GYZ, {1},, {0}\n".format(c55, self.__matid)
        self.apdl += "MPDATA, GXZ, {1},, {0}\n".format(c66, self.__matid)

        return self.__matid

    def setAreaPropByCoord(self, x, y, matId=1, createRandomCS=False):

        self.apdl += "ASEL,S,LOC,X,{0},{1}\n".format(0.95 * x, 1.05 * x)
        self.apdl += "ASEL,R,LOC,Y,{0},{1}\n".format(0.95 * y, 1.05 * y)
        csid = 0

        if createRandomCS:
            self.__csid += 1
            self.apdl += "LOCAL, {2}, 0, {0}, {1}, 0, RAND(0, 360),, , 1, 1,\n".format(x, y, self.__csid)
            csid = self.__csid

        self.apdl += "AATT,       {0}, ,   1,       {1},\n".format(matId, csid)
        self.apdl += "ASEL,S, , ,all \n"
        self.apdl += "CSYS,0  \n"

    def delOuterArea(self, x1, y1, x2, y2):
        self.apdl += "ASEL,S,LOC,X,{0},{1}\n".format(x1, x2)
        self.apdl += "ASEL,R,LOC,Y,{0},{1}\n".format(y1, y2)
        self.apdl += "ASEL, INVE\n"
        self.apdl += "ADELE,ALL, , ,1\n"
        self.apdl += "ASEL,S, , ,all\n"

    def setCirlceAreaMatProps(self, rad, matId):
        self.apdl += """
        CSYS,1  
        ASEL,S,LOC,X,0,{0}   
        CSYS,0
        AATT,       {1},
        ASEL,S, , ,all
        """.format(rad, matId)

    def setAreaProps(self, arealimit, matId=1):

        self.apdl += """
            nextarea = 0
            maxarea = {1}
            bigareaid = 0
            csid=12
            *get, numarea, area, 0, count
            
            *do, i, 1, numarea, 1
                asum,
                *get, nextarea, area, nextarea, NXTH
                *get, area_val, area, nextarea, area
                *if,area_val,gt,maxarea,then
                    bigareaid = nextarea
                    *CYCLE
                *endif
                
                ASEL,S, , , nextarea 
                asum,
                *get, cx, area, 0, cent, x
                *get, cy, area, 0, cent, y
                LOCAL, csid, 0, cx, cy, 0, RAND(0, 360),, , 1, 1,
                AATT,       {0}, ,   1, csid
                csid = csid + 1
                ASEL,S, , ,all
                CSYS, 0
            *enddo
        \n""".format(matId, arealimit)

    def mesh(self, smartsize=1):
        self.apdl += "SMRT, {0}\n".format(smartsize)
        self.apdl += "AMESH, all\n"

    def applyTensX(self, x1, y1, x2, y2, eps=0.1):

        self.prep7()
        self.clearbc()

        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x1)
        self.apdl += "DL, ALL, ,UX,0\n"

        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x2)
        self.apdl += "DL, ALL, ,UX,{0}\n".format(eps * x2)
        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y1)
        self.apdl += "LSEL,A,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UY,0\n"
        self.apdl += "LSEL,S, , ,all\n"

        self.solve()
        self.post()
        self.prep7()

        self.apdl += """
                   mat_s(1,1) = SXX0
                   mat_s(1,2) = SYY0
                   mat_s(1,3) = 0
                   vec_b(1) = EXX0\n"""

    def applyTensY(self, x1, y1, x2, y2, eps=0.1):

        self.prep7()
        self.clearbc()

        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y1)
        self.apdl += "DL, ALL, ,UY,0\n"
        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UY,{0}\n".format(eps * y2)
        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x1)
        self.apdl += "LSEL,A,LOC,X,{0}\n".format(x2)
        self.apdl += "DL, ALL, ,UX,0\n"
        self.apdl += "LSEL,S, , ,all\n"

        self.solve()
        self.post()
        self.prep7()
        self.apdl += """
            mat_s(2,1) = 0
            mat_s(2,2) = SXX0
            mat_s(2,3) = SYY0
            vec_b(2) = EYY0\n"""

    def applyTensXandY(self, x1, y1, x2, y2, eps=0.1):

        self.prep7()
        self.clearbc()

        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y1)
        self.apdl += "DL, ALL, ,UY,0\n"
        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UY,{0}\n".format(eps * y2)

        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x1)
        self.apdl += "DL, ALL, ,UX,0\n"
        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x2)
        self.apdl += "DL, ALL, ,UX,{0}\n".format(eps * x2)

        self.apdl += "LSEL,S, , ,all\n"

        self.solve()
        self.post()
        self.prep7()

        self.apdl += """
            mat_s(3,1) = 0
            mat_s(3,2) = SXX0
            mat_s(3,3) = SYY0
            vec_b(3)=EYY0\n"""

    def applyTensXandY(self, x1, y1, x2, y2, epsx, epsy):

        self.prep7()
        self.clearbc()

        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y1)
        self.apdl += "DL, ALL, ,UY,0\n"
        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UY,{0}\n".format(epsy * y2)

        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x1)
        self.apdl += "DL, ALL, ,UX,0\n"
        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x2)
        self.apdl += "DL, ALL, ,UX,{0}\n".format(epsx * x2)

        self.apdl += "LSEL,S, , ,all\n"

        self.solve()
        self.post()
        self.prep7()

        self.apdl += """
            mat_s(3,1) = 0
            mat_s(3,2) = SXX0
            mat_s(3,3) = SYY0
            vec_b(3)=EYY0\n"""

    def applyShearXY(self, x1, y1, x2, y2, eps=0.1):

        self.prep7()
        self.clearbc()

        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y1)
        self.apdl += "DL, ALL, ,ALL,0\n"
        #        self.apdl += "LSEL,S,LOC,X,{0}\n".format(x1)
        #        self.apdl += "DL, ALL, ,UY,0\n"

        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UX,{0}\n".format(eps * x2)
        self.apdl += "LSEL,S,LOC,Y,{0}\n".format(y2)
        self.apdl += "DL, ALL, ,UY,{0}\n".format(0)

        self.apdl += "LSEL,S, , ,all\n"

        self.solve()
        self.post()
        self.prep7()
        self.apdl += "GXYe = SXY0/EXY0\n"

    def solve(self):
        self.apdl += "FINISH\n"
        self.apdl += "/SOL\n"
        self.apdl += "SOLVE\n"

    def post(self):
        self.apdl += "*use,post_new1.mac\n"

    def clearbc(self):
        self.apdl += "LSCLEAR,ALL\n"

    def prep7(self):
        self.apdl += "FINISH\n"
        self.apdl += "/prep7\n"

    def eof(self):
        self.apdl += "/eof\n"

    def getAVGStressAndStrains(self):
        res = np.loadtxt(self.projdir + '/sigma_out.csv', dtype=float, delimiter=';', skiprows=1, max_rows=1)
        stress = np.zeros((3, 3))
        strain = np.zeros((3, 3))

        stress[0, 0] = res[0]
        stress[1, 1] = res[1]
        stress[2, 2] = res[2]
        stress[0, 1] = stress[1, 0] = res[3]
        stress[0, 2] = stress[2, 0] = res[4]
        stress[1, 2] = stress[2, 1] = res[5]

        strain[0, 0] = res[6]
        strain[1, 1] = res[7]
        strain[2, 2] = res[8]
        strain[0, 1] = strain[1, 0] = res[9]
        strain[0, 2] = strain[2, 0] = res[10]
        strain[1, 2] = strain[2, 1] = res[11]

        return stress, strain

    def getMaxStressForEachMaterial(self):
        res = np.loadtxt(self.projdir + '/max_stress_out.csv', dtype=float, delimiter=';', max_rows=self.__matid)
        return res

    def saveMaxStressForEachMaterial(self):
        self.apdl += """
        /post1
        *cfopen,'max_stress_out',csv,,

        *do,i,1,{0}, 1 
            ESEL,S,MAT,,i
            /SHOW, PNG 
            PLNSOL, S,EQV, 0, 1.0
            /SHOW,CLOSE
            *GET, MaxStressPar, PLNSOL, 0, MAX
            *vwrite,MaxStressPar,
%G
        *enddo
        *cfclos
        esel,all
        """.format(self.__matid)

        pass

    def precessElasticConstants(self):
        self.apdl += """
        *PRINT,mat_s
        *PRINT,vec_b

        *LSENGINE,LAPACK,MyBcsSolver,mat_s

        *LSFACTOR,MyBcsSolver
        *LSBAC,MyBcsSolver,vec_b,vec_x 

        *PRINT,vec_x

        Exe = 1.0/ vec_x(1)
        Eye = 1.0/ vec_x(3)

        nuxye = -vec_x(2)*Exe
        nuyxe = -vec_x(2)*Eye

        GXYe = SXY0/EXY0 
        
        *cfopen,'out',txt

        *vwrite, Exe, Eye, nuxye, nuyxe, GXYe 
        %/Ex=%E%/Ey=%E%/nuxy=%G%/nuyx=%G%/GXY=%E%/
        *cfclose
        \n"""

    __post_new = """

*CREATE, post_new1, mac 

FINISH
/post1

set,last
allsel,all


/SHOW,PNG,,0
PNGR,COMP,1,-1  
PNGR,ORIENT,HORIZ   
PNGR,COLOR,2
PNGR,TMOD,1 
/GFILE,1200,
!*  
/CMAP,_TEMPCMAP_,CMP,,SAVE  
/RGB,INDEX,100,100,100,0
/RGB,INDEX,0,0,0,15 

!/VIEW,1,1,1,1   
!/ANG,1  
!/AUTO,1 

/EFACET,1   
PLNSOL, S, EQV, 0, 1.0

/CMAP,_TEMPCMAP_,CMP
/DELETE,_TEMPCMAP_,CMP  
/SHOW,CLOSE 


NSEL,R,S,EQV,,, ,0  

*cfopen,'sigma_out',csv,,

*if,it_num,ne,0,then

*vwrite,'sx','sy','sz','sxy','sxz','syz','ex','ey','ez','exy','exz','eyz'
%C;%C;%C;%C;%C;%C;%C;%C;%C;%C;%C;%C

it_num = 0
*endif

ETABLE, ,VOLU, ! Get element volume
ETABLE, ,S,X ! Get element stress
ETABLE, ,S,Y
ETABLE, ,S,Z
ETABLE, ,S,XY
ETABLE, ,S,XZ
ETABLE, ,S,YZ


ETABLE, ,EPTO,X ! Get element strain
ETABLE, ,EPTO,Y
ETABLE, ,EPTO,Z
ETABLE, ,EPTO,XY
ETABLE, ,EPTO,XZ
ETABLE, ,EPTO,YZ


SMULT,SXV,VOLU,SX,1,1, ! Stress by element volume
SMULT,SYV,VOLU,SY,1,1,
SMULT,SZV,VOLU,SZ,1,1,
SMULT,SXYV,VOLU,SXY,1,1,
SMULT,SXZV,VOLU,SXZ,1,1,
SMULT,SYZV,VOLU,SYZ,1,1,

SMULT,EXV,VOLU,EPTOX,1,1, ! Strain by element volume
SMULT,EYV,VOLU,EPTOY,1,1,
SMULT,EZV,VOLU,EPTOZ,1,1,
SMULT,EXYV,VOLU,EPTOXY,1,1,
SMULT,EXZV,VOLU,EPTOXZ,1,1,
SMULT,EYZV,VOLU,EPTOYZ,1,1,

SSUM

*GET,TOTVOL,SSUM,,ITEM,VOLU ! integrate stress
*GET,TOTSX ,SSUM,,ITEM,SXV
*GET,TOTSY ,SSUM,,ITEM,SYV
*GET,TOTSZ ,SSUM,,ITEM,SZV
*GET,TOTSXY ,SSUM,,ITEM,SXYV
*GET,TOTSXZ ,SSUM,,ITEM,SXZV
*GET,TOTSYZ ,SSUM,,ITEM,SYZV

*GET,TOTVOL,SSUM,,ITEM,VOLU ! integrate strain
*GET,TOTEX ,SSUM,,ITEM,EXV
*GET,TOTEY ,SSUM,,ITEM,EYV
*GET,TOTEZ ,SSUM,,ITEM,EZV
*GET,TOTEXY ,SSUM,,ITEM,EXYV
*GET,TOTEXZ ,SSUM,,ITEM,EXZV
*GET,TOTEYZ ,SSUM,,ITEM,EYZV


SXX0 = TOTSX/TOTVOL ! compute average stress
SYY0 = TOTSY/TOTVOL
SZZ0 = TOTSZ/TOTVOL
SXY0 = TOTSXY/TOTVOL
SXZ0 = TOTSXZ/TOTVOL
SYZ0 = TOTSYZ/TOTVOL

EXX0 = TOTEX/TOTVOL ! compute average strain
EYY0 = TOTEY/TOTVOL
EZZ0 = TOTEZ/TOTVOL
EXY0 = TOTEXY/TOTVOL
EXZ0 = TOTEXZ/TOTVOL
EYZ0 = TOTEYZ/TOTVOL


*vwrite,SXX0 ,SYY0 ,SZZ0 ,SXY0 ,SXZ0 ,SYZ0, EXX0 ,EYY0 ,EZZ0 ,EXY0 ,EXZ0 ,EYZ0 
%E;%E;%E;%E;%E;%E;%E;%E;%E;%E;%E;%E

*cfclos


allsel,all

*END  
"""
