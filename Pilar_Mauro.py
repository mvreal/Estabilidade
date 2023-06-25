import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ler_arquivo(arquivo):
    
#
# Data input
#
# Data are read from a Excel datasheet called dados_pilar.xlsx
#
    path_arq = "c://users//mauro//onedrive//pilar//"
    with open(path_arq + arquivo + '.xlsx', 'rb') as target:
        sheet =  pd.read_excel(target, sheet_name='Planilha1')
        data  =  sheet.values
        
    #    f.write("Arquivo de dados: " + arquivoselecionado + "\n")
        
        ivinc, nsec, compr, fi = int(data[0,0]), int(data[0,1]), float(data[0,2]), float(data[0,3])

    #    f.write("{0:8d} {1:8d} {2:8.1f} {3:8.1f}\n".format(ivinc, nsec, compr, fi))
        
        ta, fyk, es, fck = str(data[0,4]),float(data[0,5]),float(data[0,6]),float(data[0,7])
        fyd = fyk / 1.15
        fcd = 1.1 * fck / 1.4 / 0.85

    #    f.write("{0:8s} {1:8.1f} {2:8.1f} {3:8.1f}\n".format(ta, fyk, es, fck))
        
        np0, nb0, pv, ph, mm, pp = [], [], [], [], [], []
        xp0, yp0, xb0, yb0, ass0 = [], [], [], [], []
        step = 0
        k = 0
        for i in range(nsec):
             
            np0i, nb0i = int(data[k,8]),int(data[k,9])
            step = np.max([np0i,nb0i]) 
        #    f.write("{0:8d} {1:8d}\n".format(np0i,nb0i))
            np0.append(np0i)
            nb0.append(nb0i)
            
            pvi, phi, mmi, ppi = float(data[k,10]),float(data[k,11]),float(data[k,12]),float(data[k,13])
        #    f.write("{0:8.1f} {1:8.1f} {2:8.1f} {3:8.1f}\n".format(pvi, phi, mmi, ppi))
            pv.append(pvi)
            ph.append(phi)
            mm.append(mmi)
            pp.append(ppi)
            
            xp0i, yp0i = [], []
            for j in range(np0i):
                xj, yj = float(data[k+j,14]),float(data[k+j,15])
            #    f.write("{0:8.1f} {1:8.1f}\n".format(xj,yj))
                xp0i.append(xj)
                yp0i.append(yj)
            xp0.append(xp0i)
            yp0.append(yp0i)
            
            xb0i, yb0i, ass0i = [], [], []
            for j in range(nb0i):
                xj, yj, assj = float(data[k+j,16]),float(data[k+j,17]),float(data[k+j,18])
            #    f.write("{0:8.1f} {1:8.1f} {2:8.1f}\n".format(xj,yj,assj))
                xb0i.append(xj)
                yb0i.append(yj)
                ass0i.append(assj)
            xb0.append(xb0i)
            yb0.append(yb0i)
            ass0.append(ass0i)

            k = k + step 
            
            y, b0, c0 = [0] * int(nsec), [0] * int(nsec), [0] * int(nsec)
        
    return compr,nsec,fi,ta,pv,ph,mm,pp,y,np0,nb0,xp0,yp0,xb0,yb0,ass0,es,fyd,fcd,b0,c0,ivinc,data,sheet

def estavel(compr,nsec,fi,ta,pv,ph,mm,pp,y,np0,nb0,xp0,yp0,xb0,yb0,ass0,es,fyd,fcd,b0,c0,ivinc,data,sheet,arquivo):
    
    n, v, m = [0.0] * (nsec), [0.0] * (nsec), [0.0] * (nsec)
    dx = compr / (nsec - 1)
    xs = [0.0]*(nsec)
    isec = [0]*(nsec)
    Count = 0
    while True:
        Count += 1
        n[0] = -pv[0]
        v[0] = -ph[0]
        m[0] = mm[0]
        y0 = [0] * (nsec)
        for i in range(nsec - 1):
            isec[i+1] = i + 1
            xs[i+1] = xs[i] + dx
            n[i+1] = (n[i] - pv[i + 1])
            v[i+1] = (v[i] - (pp[i] + pp[i + 1]) / 2 * dx - ph[i + 1])
            m[i+1] = (m[i] + v[i] * dx - n[i] * (y[i + 1] - y[i]) + mm[i + 1] - pp[i] * dx ** 2 / 2 - (pp[i + 1] - pp[i]) * dx ** 2 / 6)
        
        print("iteração: ",Count)

        for i in range(int(nsec)):
            y0[i] = y[i]
            np = np0[i]
            nb = nb0[i]

            xp = [xp0[i][j] for j in range(int(np))]
            yp = [yp0[i][j] for j in range(int(np))]
            xb = [xb0[i][j] for j in range(int(nb))]
            yb = [yb0[i][j] for j in range(int(nb))]
            ass = [ass0[i][j] for j in range(int(nb))]
            b0[i],c0[i] = curvatura(fi, ta, es, fyd, fcd, np, xp, yp, nb, xb, yb, ass, b0[i], c0[i], n[i], m[i])
        
        w = [0] * int(nsec)
        w[0] = dx / 12 * (3.5 * b0[0] + 3 * b0[1] - 0.5 * b0[2])
        w[-1] = dx / 12 * (3.5 * b0[-1] + 3 * b0[-2] - 0.5 * b0[-3])

        if ivinc == 1:
            r0 = w[0] + w[-1]
            m0 = -compr * w[-1]
        else:
            r0 = w[0]
            m0 = 0

        for i in range(1, int(nsec) - 1):
            w[i] = dx / 12 * (b0[i - 1] + 10 * b0[i] + b0[i + 1])
            if ivinc == 1:
                r0 += w[i]
                m0 -= dx * i * w[i]
            else:
                r0 += w[i] * (compr - i * dx) / compr
        
        y[0] = m0
        t = [0] * int(nsec)
        t[0] = r0 - w[0]
        for i in range(1, int(nsec)-1):
            y[i] = y[i - 1] + dx * t[i - 1]
            t[i] = t[i - 1] - w[i]

        a = 0
        b = 0

        for i in range(0,int(nsec)):
            a += y[i] ** 2
            b += (y[i] - y0[i]) ** 2

        if Count > 200 or b / a <= 1e-20:
            data[0,19] = Count
            data[0:nsec,20] = isec
            data[0:nsec,21] = xs
            data[0:nsec,22] = y
            data[0:nsec,23] = n
            data[0:nsec,24] = v
            data[0:nsec,25] = m
            data[0:nsec,26] = b0
        
            #
            # Saving results in a Excel datasheet: resultadosFORM.xlsx
            #

            dfres = pd.DataFrame(data,columns=sheet.columns)
            
            
            #
            path_arq = "c://users//mauro//onedrive//pilar//"
            with pd.ExcelWriter(path_arq + 'resultados_' + arquivo + '.xlsx' , engine='openpyxl') as writer:
            #    sheet.to_excel(writer, sheet_name='Planilha1', index=False)
                dfres.to_excel(writer, sheet_name='Planilha1', index=False,header='True', columns=sheet.columns)
                
            # create a figure window
            fig = plt.figure(figsize=(10, 8))
            # subplot: deslocamento y
            ax1 = fig.add_subplot(2, 2, 1)
            ax1.plot(xs,y)
            ax1.grid()
            ax1.set_xlabel('coordenada x (cm)')
            ax1.set_ylabel('deslocamento y(x) (cm)')
            ax1.set_title('Deslocamento y')
            # subplot: esforço cortante V
            ax2 = fig.add_subplot(2, 2, 2)
            ax2.plot(xs,v)
            ax2.grid()
            ax2.set_xlabel('coordenada x (cm)')
            ax2.set_ylabel('esforço cortante V(x) (kN)')
            ax2.set_title('Esforço Cortante V(x)')
            # subplot: esforço normal
            ax3 = fig.add_subplot(2, 2, 3)
            ax3.plot(xs,n)
            ax3.grid()
            ax3.set_xlabel('coordenada x (cm)')
            ax3.set_ylabel('esforço normal N(x) (kN)')
            ax3.set_title('Esforço Normal N(x)')
            # subplot: wide subplot of sinc function
            ax4 = fig.add_subplot(2, 2, 4)
            ax4.plot(xs,m)
            ax4.grid()
            ax4.set_xlabel('coordenada x (cm)')
            ax4.set_ylabel('momento fletor M(x) (kN.cm)')
            ax4.set_title('Momento Fletor M(x)')
            # plotting
            fig.tight_layout()
            fig.show()

            

            break

def curvatura(fi, ta, es, fyd, fcd, np, xp, yp, nb, xb, yb, ass, b, c, na, maxx):
    #curvatura(fi, ta, es, fyd, fcd, np, xp, yp, nb, xb, yb, ass, b0[i], c0[i], n[i], m[i])
    k = [0] * 3
    u = [0] * 2
    p = [0] * 2
    count = 0
    ruptura = False
    mrx = 0
    nr = 0

    while True:
        count += 1
        mrx, nr, k, ruptura = esforcos(b,c,fi,np,yp,nb,yb,ta,ass,xp)
        if ruptura:
            b /= 2
            c /= 2
            continue
        
        p[0] = maxx - mrx
        p[1] = na - nr
        solve(k, u, p)

        if count > 200:
            break
        
        aux = (math.sin(count * math.pi / 200.00) ** 0.1)
        b += u[0] * aux
        c += u[1] * aux

        if (p[0] ** 2 + p[1] ** 2) / (maxx ** 2 + na ** 2) < 1E-20:
            break

    if count > 200:
        print(">> RUPTURA <<")

    return b,c

def difer(i, y01, y12, xp, yp, eps0, eps1):
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
    return x1i, y1i, x2i, y2i, x1ii, y1ii, x2ii, y2ii

def esforcos(b,c,fi,np,yp,nb,yb,ta,ass,xp):
    #mrx, nr, k, ruptura = esforcos(b,c,fi,np,yp,nb,yb,ta,ass,xp)
    bb = b / (1 + fi)
    cc = c / (1 + fi)
    k = [0] * 3
    nr = 0
    mrx = 0
    ys = -1E+20
    yi = 1E+20
    ruptura = False
    comprimida = True

    for i in range(int(np)):
        
        eps = b * yp[i] + c
        if eps < -0.0035:
            ruptura = True
        if eps > 0:
            comprimida = False
        if yp[i] > ys:
            ys = yp[i]
        if yp[i] < yi:
            yi = yp[i]
    
    epss = b * ys + c
    epsi = b * yi + c

    if comprimida:
        eps = b * (ys / 1.75 + yi / 2.33333333333333) + c
        if eps < -0.002:
            ruptura = True

    for i in range(int(nb)):
        epsb = b * yb[i] + c

        if epsb > 0.01:
            ruptura = True

        sig, et = aco(epsb,ta,fyd,es)
        aux = ass[i] * et
        k[0] += aux * yb[i] * yb[i]
        k[1] += aux * yb[i]
        k[2] += aux
        aux = ass[i] * sig
        mrx += aux * yb[i]
        nr += aux

    if ruptura:
        return mrx, nr, k, ruptura
    
    if abs(b) < 0.0000000001:
        if c >= 0:
            return mrx, nr, k, ruptura
        k, mrx, nr = centra(np, fcd, bb, eps, cc, xp, yp, nr, mrx, k)
    else:
        y01 = -c / b
        y12 = (-0.002 - c) / b
        for i in range(int(np) - 1):

            eps0 = bb * yp[i] + cc
            eps1 = bb * yp[i + 1] + cc

            if abs(eps0 - eps1) > 1E-19:
                x1i, y1i, x2i, y2i, x1ii, y1ii, x2ii, y2ii = difer(i, y01, y12, xp, yp, eps0, eps1)
                k, mrx, nr = regi(fcd, bb, cc, x1i, y1i, x2i, y2i, nr, mrx, k)
                mrx, nr = regii(fcd, x1ii, y1ii, x2ii, y2ii, nr, mrx)
    return mrx, nr, k, ruptura

def centra(np, fcd, b, c, eps, xp, yp, nr, mrx, k):
    #centra(np, fcd, bb, cc, eps, xp, yp, nr, mrx, k)
    # ANTES NA CHAMADA, NO LUGAR DO EPS ELE PASSAVA A "CC" (x2)
    if eps > 0:
        return
    if eps < -0.002:
        for i in range(np-1):
            j = i + 1
            mrx, nr = regii(fcd, xp[i], yp[i], xp[j], yp[j], nr, mrx)
    else:
        for i in range(np-1):
            j = i + 1
            k, mrx, nr = regi(fcd, b, c, xp[i], yp[i], xp[j], yp[j], nr, mrx, k)
    return k, mrx, nr

def regi(fcd, b, c, x1, y1, x2, y2, nr, mrx, k):
    #regi(fcd, b, c, xp[i], yp[i], xp[j], yp[j], nr, mrx, k)
    #regi(fcd, bb, cc, x1i, y1i, x2i, y2i, nr, mrx, k)
    if x1 == 0 and y1 == 0 and x2 == 0 and y2 == 0:
        return k,mrx,nr
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
    nr = nr + sigma * (d0 * g00 + d1 * g01 + d2 * g02)
    mrx = mrx + sigma * (d0 * g01 + d1 * g02 + d2 * g03)
    k[0] = k[0] + sigma * (cle * g02 + ble * g03)
    k[1] = k[1] + sigma * (cle * g01 + ble * g02)
    k[2] = k[2] + sigma * (cle * g00 + ble * g01)
    return k, mrx, nr

def regii(fcd, x1, y1, x2, y2, nr, mrx):
    #regii(fcd, xp[i], yp[i], xp[j], yp[j], nr, mrx)
    #regii(fcd, x1ii, y1ii, x2ii, y2ii, nr, mrx)
    if x1 == 0 and y1 == 0 and x2 == 0 and y2 == 0:
        return mrx,nr
    dx = x2 - x1
    dy = y2 - y1
    g00 = (x1 + dx * 0.5) * dy
    g01 = (x1 * (y1 + dy * 0.5) + dx * (y1 * 0.5 + dy / 3)) * dy
    sigma = 0.85 * fcd
    nr -= sigma * g00
    mrx -= sigma * g01
    return mrx,nr
    
def aco(eps,ta,fyd,es):
    #aco(epsb,ta)
    epsyd = fyd / es
    sig = np.copysign(fyd, eps)
    et = 0
    if ta == "A" or ta == '"A"':
        if abs(eps) < abs(epsyd):
            sig = eps * es
            et = es
    else:
        epsyd += 0.002
        if abs(eps) < abs(0.7 * epsyd):
            sig = eps * epsyd
            et = es
        elif abs(eps) < abs(epsyd):
            a = 2.22222222222222E-02 / (fyd ** 2)
            b = 1 / es + 3.11111111111111E-02 / fyd
            c = 1.08888888888889E-02 - abs(eps)
            aux = math.sqrt(b ** 2 - 4 * a * c)
            sig = (-b + math.copysign(aux, eps)) / (2 * a)
            et = 1 / aux
    return sig, et

def solve(k, u, p):
    u[1] = (k[0] * p[1] - k[1] * p[0]) / (k[0] * k[2] - k[1] * k[1])
    u[0] = (p[0] - k[1] * u[1]) / k[0]

#
# Rodar o programa estável
#

#arquivo = str(input('Nome do arquivo de dados?'))
arquivo = str('exemplo3')
compr, nsec, fi, ta, pv, ph, mm, pp, y, np0, nb0, xp0, yp0, xb0, yb0, ass0, es, fyd, fcd, b0, c0, ivinc,data,sheet = ler_arquivo(arquivo)
estavel(compr,nsec,fi,ta,pv,ph,mm,pp,y,np0,nb0,xp0,yp0,xb0,yb0, ass0, es, fyd, fcd, b0, c0, ivinc,data,sheet,arquivo)