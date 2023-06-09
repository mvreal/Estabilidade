import math
import numpy as np

class PilarEsbelto():

    ivinc = 0
    nsec = 0
    compr = 0.00
    fck = 0.00
    fcd = 0.00
    fyk = 0.00
    fyd = 0.00
    ta = " "
    es = 0.00
    fi = 0.00
    ph = []
    pv = []
    mm = []
    pp = []
    np0 = []
    nb0 = []
    xp0 = []
    yp0 = []
    xb0 = []
    yb0 = []
    ass0 = []
    y = []
    b0 = []
    c0 = []
    y0 = []

    def __init__(self): 
        pass
        
                                                                                       
    def ler_arquivo(self, arquivoselecionado):

        with open(arquivoselecionado, "r") as entrada:
            arquivo = "saida.out"
            with open(arquivo, "w") as saida:
                saida.write("Arquivo de dados: " + arquivoselecionado + "\n")
       
                linha = entrada.readline() 
                linha = linha.split()
                self.ivinc = int(linha[0])
                self.nsec = int(linha[1])
                self.compr, self.fi = map(float, linha[2:])
                saida.write(f"{self.ivinc} {self.nsec} {self.compr} {self.fi}\n")
       
                linha = entrada.readline() 
                linha = linha.split()
                self.ta =linha[0]
                self.fyk, self.es, self.fck = map(float, linha[1:])
                self.fyd = self.fyk / 1.15
                self.fcd = 1.1 * self.fck / 1.4 / 0.85
                saida.write(f"{self.ta} {self.fyk} {self.es} {self.fck}\n")
       
                self.np0, self.nb0, self.pv, self.ph, self.mm, self.pp = [], [], [], [], [], []
                self.xp0, self.yp0, self.xb0, self.yb0, self.ass0 = [], [], [], [], []
                for i in range(1, int(self.nsec) + 1):
                    linha = entrada.readline() 
                    linha = linha.split()
                    np0i, nb0i = map(int, linha)
                    saida.write(f"{np0i} {nb0i}\n")
                    self.np0.append(np0i)
                    self.nb0.append(nb0i)
                
                    linha = entrada.readline() 
                    linha = linha.split()
                    pvi, phi, mmi, ppi = map(float, linha)
                    saida.write(f"{pvi} {phi} {mmi} {ppi}\n")
                    self.pv.append(pvi)
                    self.ph.append(phi)
                    self.mm.append(mmi)
                    self.pp.append(ppi)
                
                    xp0i, yp0i = [], []
                    for j in range(int(np0i)):
                        linha = entrada.readline() 
                        linha = linha.split()
                        xj, yj = map(float, linha)
                        saida.write(f"{xj} {yj}\n")
                        xp0i.append(xj)
                        yp0i.append(yj)
                    self.xp0.append(xp0i)
                    self.yp0.append(yp0i)
                
                    xb0i, yb0i, ass0i = [], [], []
                    for j in range(int(nb0i)):
                        linha = entrada.readline() 
                        linha = linha.split()
                        xj, yj, assj = map(float, linha)
                        saida.write(f"{xj} {yj} {assj}\n")
                        xb0i.append(xj)
                        yb0i.append(yj)
                        ass0i.append(assj)
                    self.xb0.append(xb0i)
                    self.yb0.append(yb0i)
                    self.ass0.append(ass0i)
                
        self.y, self.b0, self.c0 = np.zeros(self.nsec), np.zeros(self.nsec), np.zeros(self.nsec)
        
    def curvatura(self, np, xp, yp, nb, xb, yb, ass):
            k = [0] * 3
            u = [0] * 2
            p = [0] * 2
            count = 0

            while True:
                count += 1
                esforcos(fi, ta, es, fyd, fcd, np, xp, yp, nb, xb, yb, ass, b, c, epss, epsi, nr, mrx, k, ruptura)
                
                if ruptura == 0:
                    b /= 2
                    c /= 2
                    continue
                
                p[0] = max - mrx
                p[1] = na - nr
                solve(k, u, p)

                if count > 200:
                    break
                
                aux = (sin(count * 3.1415 / 200) ** 0.1)
                b += u[0] * aux
                c += u[1] * aux

                if (p[0] ** 2 + p[1] ** 2) / (max ** 2 + na ** 2) < 1E-20:
                    break

            if count > 200:
                print(">>>>>>>>>>>>>>> ruptura")
                # Stop (if you want to stop execution at this point, uncomment this line)


    def estavel(self):
        dx = self.compr / (self.nsec - 1)
        Count = 0
        while True:
            Count += 1
            n = [-self.pv[0]]
            v = [-self.ph[0]]
            m = [self.mm[0]]
            for i in range(self.nsec - 1):
                n.append(n[i] - self.pv[i + 1])
                v.append(v[i] - (self.pp[i] + self.pp[i + 1]) / 2 * dx - self.ph[i + 1])
                m.append(m[i] + v[i] * dx - n[i] * (self.y[i + 1] - self.y[i]) + self.mm[i + 1] - self.pp[i] * dx ** 2 / 2 - (self.pp[i + 1] - self.pp[i]) * dx ** 2 / 6)
            for i in range(self.nsec):
                self.y0[i] = self.y[i]
                np = self.np0[i]
                nb = self.nb0[i]
                xp = [self.xp0[i][j] for j in range(np)]
                yp = [self.yp0[i][j] for j in range(np)]
                xb = [self.xb0[i][j] for j in range(nb)]
                yb = [self.yb0[i][j] for j in range(nb)]
                ass = [self.ass0[i][j] for j in range(nb)]
                n[i], m[i], b, c = curvatura(self, np, xp, yp, nb, xb, yb, ass)
            w = [0] * self.nsec
            w[0] = dx / 12 * (3.5 * self.b0[0] + 3 * self.b0[1] - 0.5 * self.b0[2])
            w[-1] = dx / 12 * (3.5 * self.b0[-1] + 3 * self.b0[-2] - 0.5 * self.b0[-3])
            if self.ivinc == 1:
                r0 = w[0] + w[-1]
                m0 = -self.compr * w[-1]
            else:
                r0 = w[0]
                m0 = 0
            for i in range(1, nsec - 1):
                w[i] = dx / 12 * (self.b0[i - 1] + 10 * self.b0[i] + self.b0[i + 1])
                if self.ivinc == 1:
                    r0 += w[i]
                    m0 -= dx * (i - 1) * w[i]
                else:
                    r0 += w[i] * (self.compr - (i - 1) * dx) / self.compr
            y[0] = m0
            t = [0] * self.nsec
            t[0] = r0 - w[0]
            for i in range(1, self.nsec - 1):
                self.y[i] = self.y[i - 1] + dx * t[i - 1]
                t[i] = t[i - 1] - w[i]
            A = 0
            b = 0
            for i in range(self.nsec - 1):
                A += self.y[i] ** 2
                b += (self.y[i] - self.y0[i]) ** 2
            if Count > 200 or b / A <= 1e-20:
                break
        


    


    def difer(i, y01, y12, xp, yp, eps0, eps1, x1i, y1i, x2i, y2i, x1ii, y1ii, x2ii, y2ii):
        t01 = 0
        t12 = 0
        x1i = 0
        y1i = 0
        x2i = 0
        y2i = 0
        x1ii = 0
        y1ii = 0
        x2ii = 0
        y2ii = 0
        
        i2 = i + 1
        dy = yp[i2] - yp[i]
        dxdy = (xp[i2] - xp[i]) / dy
        dum1 = y01 - yp[i]
        dum2 = y12 - yp[i]
        x01 = xp[i] + dum1 * dxdy
        x12 = xp[i] + dum2 * dxdy
        dy01 = dum1 / dy
        dy12 = dum2 / dy
        
        if dy01 > 0 and dy01 < 1:
            t01 = 1
        if dy12 > 0 and dy12 < 1:
            t12 = 1
        
        if eps0 < eps1:
            t01 = -t01
            t12 = -t12
            
        if t01 == 0 and t12 == 0:
            if eps0 < 0:
                if eps0 > -0.002:
                    x1i = xp[i]
                    y1i = yp[i]
                    x2i = xp[i2]
                    y2i = yp[i2]
                else:
                    x1ii = xp[i]
                    y1ii = yp[i]
                    x2ii = xp[i2]
                    y2ii = yp[i2]
        
        else:
            if t01 == 1:
                x1i = x01
                y1i = y01
                if t12 == 1:
                    x2i = x12
                    y2i = y12
                    x1ii = x12
                    y1ii = y12
                    x2ii = xp[i2]
                    y2ii = yp[i2]
                else:
                    x2i = xp[i2]
                    y2i = yp[i2]
            
            else:
                if t01 == -1:
                    x2i = x01
                    y2i = y01
                    if t12 == -1:
                        x1i = x12
                        y1i = y12
                        x2ii = x12
                        y2ii = y12
                        x1ii = xp[i]
                        y1ii = yp[i]
                    else:
                        x1i = xp[i]
                        y1i = yp[i]
                
                else:
                    if t12 == 1:
                        x1i = xp[i]
                        y1i = yp[i]
                        x2i = x12
                        y2i = y12
                        x1ii = x12
                        y1ii = y12
                        x2ii = xp[i2]
                        y2ii = yp[i2]
                    else:
                        x1i = x12
                        y1i = y12
                        x2i = xp[i2]
                        y2i = yp[i2]
                        x1ii = xp[i]
                        y1ii = yp[i]
                        x2ii = x12

    def esforcos(fi, ta, es, fyd, fcd, np, xp, yp, nb, xb, yb, ass, b, c, epss, epsi, nr, mrx, k, ruptura):
        bb = b / (1 + fi)
        cc = c / (1 + fi)
        sim = 0
        nao = 1
        k[0], k[1], k[2] = 0, 0, 0
        nr = 0
        mrx = 0
        ys = -1E+20
        yi = 1E+20
        ruptura = nao
        comprimida = sim
        for i in range(np):
            eps = b * yp[i] + c
            if eps < -0.0035:
                ruptura = sim
            if eps > 0:
                comprimida = nao
            if yp[i] > ys:
                ys = yp[i]
            if yp[i] < yi:
                yi = yp[i]
        epss = b * ys + c
        epsi = b * yi + c
        if comprimida == sim:
            eps = b * (ys / 1.75 + yi / 2.33333333333333) + c
            if eps < -0.002:
                ruptura = sim
        for i in range(nb):
            epsb = b * yb[i] + c
            if epsb > 0.01:
                ruptura = sim
            aco(epsb, ta, fyd, es, et, sig)
            aux = ass[i] * et
            k[0] += aux * yb[i] * yb[i]
            k[1] += aux * yb[i]
            k[2] += aux
            aux = ass[i] * sig
            mrx += aux * yb[i]
            nr += aux

        if ruptura == sim:
            return
        if abs(b) < 0.0000000001:
            if c >= 0:
                return
            centra(np, fcd, bb, cc, cc, xp, yp, nr, mrx, k)
        else:
            y01 = -c / b
            y12 = (-0.002 - c) / b
            for i in range(np - 1):
                eps0 = bb * yp[i] + cc
                eps1 = bb * yp[i + 1] + cc
                if abs(eps0 - eps1) > 1E-19:
                    difer(i, y01, y12, xp, yp, eps0, eps1, x1i, y1i, x2i, y2i, x1ii, y1ii, x2ii, y2ii)
                    regi(fcd, bb, cc, x1i, y1i, x2i, y2i, nr, mrx, k)
                    regii(fcd, x1ii, y1ii, x2ii, y2ii, nr, mrx)

    def centra(np, fcd, b, c, eps, xp, yp, nr, mrx, k):
        if eps > 0:
            return
        if eps < -0.002:
            for i in range(1, np):
                j = i + 1
                regi(fcd, b, c, xp[i-1], yp[i-1], xp[j-1], yp[j-1], nr, mrx, k)
        else:
            for i in range(1, np):
                j = i + 1
                regii(fcd, xp[i-1], yp[i-1], xp[j-1], yp[j-1], nr, mrx)

    def regi(fcd, b, c, x1, y1, x2, y2, nr, mrx, k):
        if x1 == 0 and y1 == 0 and x2 == 0 and y2 == 0:
            return
        cle = 500 * c + 1
        ble = 500 * b
        d0 = c * (1 + 250 * c)
        d1 = b * cle
        d2 = b * ble * 0.5
        sigma = 850 * fcd
        dx = x2 - x1
        dy = y2 - y1
        dy1 = dy * 0.5
        dy2 = dy * dy
        dy3 = dy2 * dy
        dx2 = dx * dx
        g00 = (x1 + dx * 0.5) * dy
        g01 = (x1 * (y1 + dy1) + dx * (y1 * 0.5 + dy / 3)) * dy
        g02 = (x1 * (y1 * (dy + y1) + dy2 / 3) + dx * (y1 * (y1 * 0.5 + dy / 1.5) + dy2 * 0.25)) * dy
        g03 = (x1 * (y1 * (dy2 + y1 * (1.5 * dy + y1)) + dy3 * 0.25) + dx * (y1 * (0.75 * dy2 + y1 * (dy + y1 * 0.5)) + dy3 * 0.2)) * dy
        nr[0] = nr[0] + sigma * (d0 * g00 + d1 * g01 + d2 * g02)
        mrx[0] = mrx[0] + sigma * (d0 * g01 + d1 * g02 + d2 * g03)
        k[0] = k[0] + sigma * (cle * g02 + ble * g03)
        k[1] = k[1] + sigma * (cle * g01 + ble * g02)
        k[2] = k[2] + sigma * (cle * g00 + ble * g01)


    def regii(fcd, x1, y1, x2, y2, nr, mrx):
        if x1 == 0 and y1 == 0 and x2 == 0 and y2 == 0:
            return
        
        dx = x2 - x1
        dy = y2 - y1
        g00 = (x1 + dx * 0.5) * dy
        g01 = (x1 * (y1 + dy * 0.5) + dx * (y1 * 0.5 + dy / 3)) * dy
        sigma = 0.85 * fcd
        nr -= sigma * g00
        mrx -= sigma * g01
        
        
    def aco(eps, ta, fyd, es):
        epsyd = fyd / es
        sig = math.copysign(fyd, eps)
        et = 0
        
        if ta == "A":
            if abs(eps) < epsyd:
                sig = eps * es
                et = es
        else:
            epsyd += 0.002
            if abs(eps) < 0.7 * epsyd:
                sig = eps * epsyd
                et = es
            elif abs(eps) < epsyd:
                A = 2.22222222222222E-02 / (fyd ** 2)
                b = 1 / es + 3.11111111111111E-02 / fyd
                c = 1.08888888888889E-02 - abs(eps)
                aux = math.sqrt(b ** 2 - 4 * A * c)
                sig = (-b + math.copysign(aux, eps)) / (2 * A)
                et = 1 / aux
        
        return sig, et    

    def solve(k, u, p):
        
        u[2] = (k[1] * p[2] - k[2] * p[1]) / (k[1] * k[3] - k[2] * k[2])
        u[1] = (p[1] - k[2] * u[2]) / k[1]

