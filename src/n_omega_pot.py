import sys, re, os, glob
import numpy as np
import pickle
import struct


binSize = 31
ainv = 2.333 # GeV
# full stat
M_N = 0.95469/ainv # in lat. unit
M_OMEGA = 1.71153/ainv # in lat. unit

M_red =  M_N * M_OMEGA/(M_N + M_OMEGA)
D_fac = (M_N - M_OMEGA)/(M_N + M_OMEGA)
Ns = 96

Header = '/nas/nasbee03/iritani/kei/results.t_mod02_001.rot_tot.ave.full.stat/'
CorrRdir = '/nas/nasbee03/iritani/kei/results.t_mod02_001.rot_tot.ave.full.stat/'

def main():
    for it in [10, 11, 12, 13, 14, 15]:
        nbs = NBS_pot_DO_Proj(it, binSize, spin=2, reset=False)

        nbs.calc_pot()

class NBS_pot_DO_Proj(object):
    """ potential for DO """
    def __init__(self, it, binSize, spin=2, confMax=None, reset=False): 
        self.it = it
        self.binSize = binSize
        self.spin = spin
        self.confMax = confMax

        self.C_N = Corr_2pt_Baryon('proton_CG05_CG05', binSize, Header, CorrRdir, confMax=confMax)
        self.C_Omega = Corr_2pt_Baryon('OMEGA', binSize, Header, CorrRdir, confMax=confMax)

        self.nconf = self.C_N.nconf
        self.binNum = self.C_N.binNum


        self.Rcorr_tm_jk = self.loadRcorr(it-1, reset)
        self.Rcorr_t_jk  = self.loadRcorr(it, reset)
        self.Rcorr_tp_jk = self.loadRcorr(it+1, reset)

    def loadRcorr(self, mit, reset=False):
        dir_name = 'results.pot.bin/spin{:d}_ave'.format(self.spin)
        if not os.path.isdir(dir_name): os.makedirs(dir_name)
    
        f_name = '{}/Rcorr_nomega_spin{:d}_t{:03d}_{:03d}conf_{:03d}bin_a1proj.bin'.format(
                    dir_name, self.spin, mit, self.binNum * self.binSize, self.binSize)

        fsize = Ns**3*self.binNum * 2

        if os.path.isfile(f_name) and reset == False:
            print('load ', f_name)
            with open(f_name, 'rb') as infile: 
                _Rcorr_jk = np.array(struct.unpack('{:d}d'.format(fsize), infile.read(8*fsize))).reshape(self.binNum,Ns,Ns,Ns,2)
                Rcorr_jk = _Rcorr_jk[:,:,:,:,0] + _Rcorr_jk[:,:,:,:,1]*1j
        else:
            print('calc ', f_name)

            wave_t_jk  = self.iowave(mit)
            Rcorr_jk = np.array([wave_t_jk[ibin,:,:,:]/self.C_N.corr_jk[ibin,mit]/self.C_Omega.corr_jk[ibin,mit]
                                for ibin in range(self.binNum)])
            with open(f_name, 'wb') as fout: 
                fout.write(bytearray(Rcorr_jk.flatten()))

        return np.real(Rcorr_jk)

    def calc_pot(self, fitr=(0,20)):
        lap = lambda vec: - 6.0*vec + (  np.roll(vec,+1,0) + np.roll(vec,-1,0) 
                                       + np.roll(vec,+1,1) + np.roll(vec,-1,1)
                                       + np.roll(vec,+1,2) + np.roll(vec,-1,2))


        print('# calc pot')
        vlap_jk = np.array([ lap(self.Rcorr_t_jk[ibin,:,:,:])/self.Rcorr_t_jk[ibin,:,:,:]/(2.0*M_red)
            for ibin in range(self.binNum)])
        vdt_jk = np.array([
            - (self.Rcorr_tp_jk[ibin,:,:,:] - self.Rcorr_tm_jk[ibin,:,:,:])/(2.0*self.Rcorr_t_jk[ibin,:,:,:])
            for ibin in range(self.binNum)])

        vdt2_jk = np.array([
              (self.Rcorr_tp_jk[ibin,:,:,:] - 2.0 * self.Rcorr_t_jk[ibin,:,:,:]
               + self.Rcorr_tm_jk[ibin,:,:,:])/(self.Rcorr_t_jk[ibin,:,:,:])*(1.0+3.0*D_fac**2)/(8.0*M_red)
            for ibin in range(self.binNum)])

        rs = np.array([np.sqrt(x**2 + y**2 + z**2) for z in range(-Ns//2,Ns//2)
                                                   for y in range(-Ns//2,Ns//2)
                                                   for x in range(-Ns//2,Ns//2)]).reshape(Ns,Ns,Ns)
        rs2 = np.roll(np.roll(np.roll(rs,Ns//2,0),Ns//2,1),Ns//2,2).flatten()
        vlap_jk = vlap_jk.reshape(self.binNum,Ns**3)
        vdt_jk = vdt_jk.reshape(self.binNum,Ns**3)
        vdt2_jk = vdt2_jk.reshape(self.binNum,Ns**3)

        vc_jk = vlap_jk + vdt_jk + vdt2_jk

        print('# save pot')
        if not os.path.isdir('results.pot.bin/spin{:d}_ave'.format(self.spin)): 
            os.makedirs('results.pot.bin/spin{:d}_ave'.format(self.spin))
        with open('results.pot.bin/spin{:d}_ave/pot_nomega_spin{:d}_t{:03d}_{:03d}conf_{:03d}bin_a1proj.bin'.format(
                        self.spin, self.spin, self.it, self.binNum * self.binSize, self.binSize), 'wb') as fout:
            fout.write(bytearray(vc_jk.flatten()))

# output potential

        uniq_a1 = [ix + Ns*(iy + Ns*iz) for iz in range(0,Ns//2+1)
                        for iy in range(iz,Ns//2+1) for ix in range(iy,Ns//2+1)]

        if not os.path.isdir('results.pot/spin{:d}_ave'.format(self.spin)): 
            os.makedirs('results.pot/spin{:d}_ave'.format(self.spin))

        with open('results.pot/spin{:d}_ave/pot_nomega_spin{:d}_t{:03d}_{:03d}conf_{:03d}bin_a1proj.dat'.format(
             self.spin, self.spin, self.it, self.binNum * self.binSize, self.binSize), 'w') as fout:

            print('# output results')

            fout.write('# r   H0  d/dt  d2/dt2 tot.\n')
            for r in uniq_a1:
                vc_av, vc_err = jack_sample_mean(vc_jk[:,r])
                vlap_av, vlap_err = jack_sample_mean(vlap_jk[:,r])
                vdt_av, vdt_err = jack_sample_mean(vdt_jk[:,r])
                vdt2_av, vdt2_err = jack_sample_mean(vdt2_jk[:,r])
                fout.write('{:e}   {:e} {:e}   {:e} {:e}   {:e} {:e}    {:e} {:e}\n'.format(rs2[r], 
                    vlap_av, vlap_err, vdt_av, vdt_err, 
                    vdt2_av, vdt2_err, vc_av, vc_err))


    def iowave(self, it):
        dir_list = os.listdir(Header + '/Proj.DOwave.dir.S3.35/spin{:d}_ave/'.format(self.spin))
        dir_list.sort()
        waves = []
        waves_bin = []
        for dir_name in dir_list[:self.confMax]:
            files = glob.glob(Header + '/Proj.DOwave.dir.S3.35/spin{:d}_ave/'.format(self.spin)
                            + dir_name + '/NBSwave.+{:03d}*'.format(it))
            for fname in files:
                with open(fname, 'rb') as infile:
                    tmpw = np.array(struct.unpack('>{:d}d'.format(Ns**3*2), infile.read(8*Ns**3*2)))
                    wave_a1proj = self.A1_projection(
                            tmpw.reshape(Ns,Ns,Ns,2)[:,:,:,0]
                            + 1j*tmpw.reshape(Ns,Ns,Ns,2)[:,:,:,1])
                    waves.append(wave_a1proj)

                if len(waves) == self.binSize:
                    print('--binning--NBSwave-- t = {:d}'.format(it))
                    waves_bin.append(np.mean(np.array(waves)[:,:,:,:], axis=0))
                    waves = []

        waves_bin = np.array(waves_bin)
        self.binNum = waves_bin.shape[0]
        wave_jk = (np.sum(waves_bin[:,:,:,:], axis=0)
                - waves_bin)/float(self.binNum - 1)
        return wave_jk

    def A1_projection(self, wave):
        wave_tmp1 = (wave[:,:,:] + np.roll(wave,-1,0)[::-1,:,:]
                    + np.roll(wave,-1,1)[:,::-1,:]
                    + np.roll(wave,-1,2)[:,:,::-1]
                    + np.roll(np.roll(wave,-1,0),-1,1)[::-1,::-1,:]
                    + np.roll(np.roll(wave,-1,1),-1,2)[:,::-1,::-1]
                    + np.roll(np.roll(wave,-1,2),-1,0)[::-1,:,::-1]
                    + np.roll(np.roll(np.roll(wave,-1,0),-1,1),-1,2)[::-1,::-1,::-1])/8.0
        wave_tmp2 = (wave_tmp1 
                    + np.swapaxes(wave_tmp1,0,1)
                    + np.swapaxes(wave_tmp1,1,2)
                    + np.swapaxes(wave_tmp1,2,0)
                    + np.swapaxes(np.swapaxes(wave_tmp1,0,1),1,2)
                    + np.swapaxes(np.swapaxes(wave_tmp1,0,2),2,1))/6.0e0

        return wave_tmp2



class Corr_2pt_Baryon(object):
    """ baryon 2pt correlator """
    def __init__(self, corr_type, binSize, header, rdir, confMax = None):
        self.corr_type = corr_type
        self.binSize = binSize
        self.rdir = rdir

        dir_list = os.listdir(header + '/correlator.PS.dir/')
        dir_list.sort()
        corrs = []
        corrs_bin = {}
        print('##Read {}'.format(corr_type))
        for fw_bw in ['', 'anti']:
            corrs = []
            corrs_bin[fw_bw] = []
            for dir_name in dir_list[:confMax]:
                files = glob.glob(header + '/correlator.PS.dir/' + dir_name + '/' 
                              + fw_bw + corr_type + '_*')
                for file in files:
#                    print('read >> ...', file[-70:])
                    corrs.append(np.loadtxt(file))
                
                    if len(corrs) == binSize:
#                    print('--binning--')
                        corrs_bin[fw_bw].append(np.mean(np.array(corrs)[:,:,1], axis=0))
                        corrs = []

        corrs_fw_bin = np.array(corrs_bin[''])
        corrs_bw_bin = np.array(corrs_bin['anti'])

        self.binNum, self.nt = corrs_fw_bin.shape
        self.nconf = self.binNum * binSize

        self.proj(corrs_fw_bin, corrs_bw_bin)


        nega_posi = ['pos', 'neg']
        ip = 0
        self.corr_jk = (np.sum(self.corrs_2pt_fb[ip][:,:], axis=0) 
                    - self.corrs_2pt_fb[ip])/float(self.binNum - 1)
        self.corr_jk_av = np.array([jack_sample_mean(self.corr_jk[:,it]) 
                                for it in np.arange(0, self.nt)])


        
    def proj(self, corrs_fw_bin, corrs_bw_bin):
        corrs_2pt_baryon = [corrs_fw_bin, corrs_fw_bin, corrs_bw_bin, corrs_bw_bin]
        
        corrs_2pt_pm_fw = self.calc_proj_pm_2pt(corrs_2pt_baryon, (0,1))
        corrs_2pt_pm_bw = self.calc_proj_pm_2pt(corrs_2pt_baryon, (1,0))
        self.corrs_2pt_fb = self.folding(corrs_2pt_pm_fw, corrs_2pt_pm_bw)

        
    def folding(self, corrs_2pt_pm_fw, corrs_2pt_pm_bw):
        SrcT = 1 # 
        Iflg_bc = 'p'
        corrs_2pt_fb = {}
        its = np.arange(self.nt)
        itrevs = (self.nt - its + 2*(SrcT - 1)) % self.nt

        for ip, sign in enumerate([+1.0, -1.0]):
            if Iflg_bc == 'd':
                corrs_2pt_fb[ip] = corrs_2pt_pm_fw[ip][:,its]
            elif Iflg_bc == 'p':
                corrs_2pt_fb[ip] = sign * 0.5 * (corrs_2pt_pm_fw[ip][:,its] 
                          + corrs_2pt_pm_bw[ip][:,itrevs])
            elif Iflg_bc == 'a':
                corrs_2pt_fb[ip] = sign * 0.5 * (corrs_2pt_pm_fw[ip][:,its]
                          - corrs_2pt_pm_bw[ip][:,itrevs])
                    
        return corrs_2pt_fb
        
    def calc_proj_pm_2pt(self, corr_2pt, mode):
# projection  ---      particle (1 + gamma_4)/2, (1 - gamma_4)/2
#             --- anti-partycle (1 - gamma_4)/2, (1 + gamma_4)/2
        corrs_2pt_parity = {}
        corrs_2pt_parity[mode[0]] = corr_2pt[0] + corr_2pt[1]
        corrs_2pt_parity[mode[1]] = corr_2pt[2] + corr_2pt[3]
        return corrs_2pt_parity



def jack_sample_mean(jack_array):
    return [np.mean(jack_array), np.sqrt(len(jack_array) - 1) * np.std(jack_array)]


if __name__ == '__main__':
    main()
