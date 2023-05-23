#merging of signal regions with given weights

from os import path
import numpy as np
from general import alp_constants as c
from alp import alp_setup as setup

class MergeInput:

	def __init__(self,exp,regions,channels_decay,channels_production):
		self._outDir = "/../../tab_decay/"
		# self._weights = [1. , 0.734375 , 1. , 0.687462687 , 0.661781609 , 1. , 1. , 1.]
		self._weightsForFiles = {}

		for chan_dec in channels_decay:
			for chan_prod in channels_production:
				files = []
				
				for reg in regions:
					file = setup.experiments[exp]+"/"+exp+'_'+chan_prod+'_'+chan_dec+reg+'.dat'
					self._outDir = path.dirname(path.realpath(__file__))+"/../../tab_decay/"
					filename_dat = self._outDir + file
					if path.exists(filename_dat):
						files.append(file)
					else:
						print('[Warning:] \t',filename_dat,'not found. Skipping..')

				if not len(files) == 0:
					self.getWeights(files)
					self.mergeRegions()

	def getWeights(self,files):
		print('[Info:] \t','Re-weighting yields for given signal regions')

		### uncomment for using fixed weights
		# weights = [1. , 0.734375 , 1. , 0.687462687 , 0.661781609 , 1. , 1. , 1.]
		# self._weightsForFiles = dict(zip(files, weights))
		#####################################

		### read weights only once, then reuse
		if not self._weightsForFiles:
			for file in files:
				self._weightsForFiles[file] = float(input(" - Enter weight for file " + file + ": "))
		else:
			if len(files) == len(list(self._weightsForFiles.values())):
				### use old values for new files
				self._weightsForFiles = dict(zip(files, list(self._weightsForFiles.values())))
				for file in files:
					print(' - Using weight for file', file, ': ',self._weightsForFiles[file])
			else:
				print('[Error:] some weight or filename not read correctly, could not merge. Exiting..')
				exit(1)

		return

	def mergeRegions(self):
		print('[Info:] \t','Merging signal regions..')
		rescaled_data_dat = None
		for file in self._weightsForFiles.keys():
			filename = self._outDir + file
			filename_dat = path.dirname(filename)

			if path.exists(filename_dat):
				
				rescaled_data_dat_in = np.loadtxt(filename)
				if rescaled_data_dat is None:
					rescaled_data_dat = self._weightsForFiles[file]*rescaled_data_dat_in
				else:
					rescaled_data_dat[:,2] += self._weightsForFiles[file]*rescaled_data_dat_in[:,2]

			else:
				print('[Error:] \t',filename_dat,'not found, could not merge. Exiting..')
				exit(1)
		
		if list(self._weightsForFiles.keys())[0].endswith('_reg1.dat'):
			outfileName = list(self._weightsForFiles.keys())[0].replace('_reg1.dat','.dat')
		else:
			print('[Error:] \t',list(self._weightsForFiles.keys())[0],'name format inconsistent. Exiting..')
			exit(1)

		np.savetxt(self._outDir + outfileName,rescaled_data_dat)
		print('[Info:] \t', 'File ' + outfileName + ' saved to ' + self._outDir)