import cPickle, base64
try:
	from SimpleSession.versions.v65 import beginRestore,\
	    registerAfterModelsCB, reportRestoreError, checkVersion
except ImportError:
	from chimera import UserError
	raise UserError('Cannot open session that was saved in a'
	    ' newer version of Chimera; update your version')
checkVersion([1, 13, 1, 41965])
import chimera
from chimera import replyobj
replyobj.status('Restoring session...', \
    blankAfter=0)
replyobj.status('Beginning session restore...', \
    blankAfter=0, secondary=True)
beginRestore()

def restoreCoreModels():
	from SimpleSession.versions.v65 import init, restoreViewer, \
	     restoreMolecules, restoreColors, restoreSurfaces, \
	     restoreVRML, restorePseudoBondGroups, restoreModelAssociations
	molInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVRFyaWJib25JbnNpZGVDb2xvcnECSwZOfYdVCWJhbGxTY2FsZXEDSwZHP9AAAAAAAAB9h1UJcG9pbnRTaXplcQRLBkc/8AAAAAAAAH2HVQVjb2xvcnEFSwZLAH1xBihLAV1xB0sBYUsCXXEISwJhSwNdcQlLA2FLBF1xCksEYUsFXXELSwVhdYdVCnJpYmJvblR5cGVxDEsGSwB9h1UKc3RpY2tTY2FsZXENSwZHP/AAAAAAAAB9h1UMbW1DSUZIZWFkZXJzcQ5dcQ8oTk5OTk5OZVUMYXJvbWF0aWNNb2RlcRBLBksBfYdVCnZkd0RlbnNpdHlxEUsGR0AUAAAAAAAAfYdVBmhpZGRlbnESSwaJfYdVDWFyb21hdGljQ29sb3JxE0sGTn2HVQ9yaWJib25TbW9vdGhpbmdxFEsGSwB9h1UJYXV0b2NoYWlucRVLBoh9h1UKcGRiVmVyc2lvbnEWSwZLAn2HVQhvcHRpb25hbHEXfXEYVQhvcGVuZWRBc3EZiIlLBihVTi9ob21lL3lhbmdqYy9naXQvY2FuLWFpLWRvL3BkYmJpbmQvYXRvbV9tYXBwaW5nL3BkYi9hdG9tLjJxYnJfbGlnYW5kLmFsaWduLnBkYk5VA1BEQnEaSwF0cRt9cRwoKFVOL2hvbWUveWFuZ2pjL2dpdC9jYW4tYWktZG8vcGRiYmluZC9hdG9tX21hcHBpbmcvcGRiL2F0b20uMmhiMV9saWdhbmQuYWxpZ24ucGRiTmgaSwF0cR1dcR4oSwBLA2UoVU4vaG9tZS95YW5namMvZ2l0L2Nhbi1haS1kby9wZGJiaW5kL2F0b21fbWFwcGluZy9wZGIvYXRvbS4ycWJwX2xpZ2FuZC5hbGlnbi5wZGJOaBpLAXRxH11xIChLAksFZXWHh3NVD2xvd2VyQ2FzZUNoYWluc3EhSwaJfYdVCWxpbmVXaWR0aHEiSwZHP/AAAAAAAAB9h1UPcmVzaWR1ZUxhYmVsUG9zcSNLBksAfYdVBG5hbWVxJEsGWBoAAABhdG9tLjJxYnBfbGlnYW5kLmFsaWduLnBkYn1xJShYGgAAAGF0b20uMmhiMV9saWdhbmQuYWxpZ24ucGRiXXEmKEsASwNlWBoAAABhdG9tLjJxYnJfbGlnYW5kLmFsaWduLnBkYl1xJyhLAUsEZXWHVQ9hcm9tYXRpY0Rpc3BsYXlxKEsGiX2HVQ9yaWJib25TdGlmZm5lc3NxKUsGRz/pmZmZmZmafYdVCnBkYkhlYWRlcnNxKl1xKyh9cSx9cS19cS5YBQAAAEhFTElYXXEvWEwAAABIRUxJWCAgICAxICAgMSBHTE4gICAgIDIxICBBUkcgICAgIDI0ICAxICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAzcTBhc31xMX1xMn1xM1gFAAAASEVMSVhdcTRYTAAAAEhFTElYICAgIDEgICAxIEdMTiAgICAgMjEgIEFSRyAgICAgMjQgIDEgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIDNxNWFzZVUDaWRzcTZLBksASwCGfXE3KEtmSwCGXXE4SwVhSwJLAIZdcTlLAmFLZUsAhl1xOksEYUsBSwCGXXE7SwFhS2RLAIZdcTxLA2F1h1UOc3VyZmFjZU9wYWNpdHlxPUsGR7/wAAAAAAAAfYdVEGFyb21hdGljTGluZVR5cGVxPksGSwJ9h1UUcmliYm9uSGlkZXNNYWluY2hhaW5xP0sGiH2HVQdkaXNwbGF5cUBLBoh9cUGJXXFCKEsBSwRlc4d1Lg=='))
	resInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQZpbnNlcnRxAksGVQEgfYdVC2ZpbGxEaXNwbGF5cQNLBol9h1UEbmFtZXEESwZYAwAAAE1PTH2HVQVjaGFpbnEFSwZYAQAAAEF9h1UOcmliYm9uRHJhd01vZGVxBksGSwJ9h1UCc3NxB0sGiYmGfYdVCG1vbGVjdWxlcQhLBksAfXEJKEsBTl1xCksBSwGGcQthhksCTl1xDEsCSwGGcQ1hhksDTl1xDksDSwGGcQ9hhksETl1xEEsESwGGcRFhhksFTl1xEksFSwGGcRNhhnWHVQtyaWJib25Db2xvcnEUSwZOfXEVKEsITl1xFksFSwGGcRdhhksGTl1xGEsDSwGGcRlhhksHTl1xGksESwGGcRthhnWHVQVsYWJlbHEcSwZYAAAAAH2HVQpsYWJlbENvbG9ycR1LBk59cR4oSwlOXXEfSwNLAYZxIGGGSwtOXXEhSwVLAYZxImGGSwpOXXEjSwRLAYZxJGGGdYdVCGZpbGxNb2RlcSVLBksBfYdVBWlzSGV0cSZLBoh9h1ULbGFiZWxPZmZzZXRxJ0sGTn2HVQhwb3NpdGlvbnEoXXEpKEsBSwGGcSpLAUsBhnErSwFLAYZxLEsBSwGGcS1LAUsBhnEuSwFLAYZxL2VVDXJpYmJvbkRpc3BsYXlxMEsGiX2HVQhvcHRpb25hbHExfVUEc3NJZHEySwZK/////32HdS4='))
	atomInfo = cPickle.loads(base64.b64decode('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'))
	bondInfo = cPickle.loads(base64.b64decode('gAJ9cQEoVQVjb2xvcnECS3xOfYdVBWF0b21zcQNdcQQoXXEFKEsZSw1lXXEGKEsMSw1lXXEHKEsMSw5lXXEIKEsMSxplXXEJKEsVSxRlXXEKKEsVSxZlXXELKEsVSxhlXXEMKEsNSxBlXXENKEsPSw5lXXEOKEsPSxBlXXEPKEsPSxFlXXEQKEsQSxNlXXERKEsRSxJlXXESKEsRSxdlXXETKEsUSxNlXXEUKEsUSxtlXXEVKEsUSxxlXXEWKEsrSyBlXXEXKEsmSyVlXXEYKEsmSydlXXEZKEsmSzllXXEaKEsnSyhlXXEbKEsnSzFlXXEcKEsoSyllXXEdKEsoSzplXXEeKEspSyplXXEfKEspSztlXXEgKEsqSyVlXXEhKEsqSzxlXXEiKEseSx1lXXEjKEseSx9lXXEkKEseSyJlXXElKEstSz1lXXEmKEstSz5lXXEnKEstSyxlXXEoKEstSy5lXXEpKEsuSy9lXXEqKEsuSzBlXXErKEsfSyBlXXEsKEsfSyxlXXEtKEsyS0BlXXEuKEsyS0FlXXEvKEsySzFlXXEwKEsySzNlXXExKEszSzRlXXEyKEszSzhlXXEzKEs0S0JlXXE0KEs0SzVlXXE1KEs1S0NlXXE2KEs1SzZlXXE3KEs2S0RlXXE4KEs2SzdlXXE5KEs3S0VlXXE6KEs3SzhlXXE7KEs4S0ZlXXE8KEsgSyFlXXE9KEshSx1lXXE+KEshSyVlXXE/KEsiSyNlXXFAKEsiSyRlXXFBKEs/SzFlXXFCKEtVS0plXXFDKEtQS09lXXFEKEtQS1FlXXFFKEtQS2xlXXFGKEtRS1JlXXFHKEtRS1tlXXFIKEtSS1NlXXFJKEtSS21lXXFKKEtTS1RlXXFLKEtTS25lXXFMKEtUS09lXXFNKEtUS29lXXFOKEtIS0dlXXFPKEtIS0llXXFQKEtIS0xlXXFRKEtXS1ZlXXFSKEtXS1hlXXFTKEtXS3BlXXFUKEtXS3FlXXFVKEtYS1llXXFWKEtYS1plXXFXKEtJS0plXXFYKEtJS1ZlXXFZKEtcS1tlXXFaKEtcS11lXXFbKEtcS2FlXXFcKEtcS3NlXXFdKEtdS15lXXFeKEtdS3RlXXFfKEtdS3VlXXFgKEteS19lXXFhKEteS3ZlXXFiKEteS3dlXXFjKEtgS19lXXFkKEtgS2FlXXFlKEtgS3hlXXFmKEtgS3llXXFnKEthS3plXXFoKEthS3tlXXFpKEtKS0tlXXFqKEtlS2JlXXFrKEtlS2ZlXXFsKEtlS3xlXXFtKEtlS31lXXFuKEtLS0dlXXFvKEtLS09lXXFwKEtmS2dlXXFxKEtmS2tlXXFyKEtnS2hlXXFzKEtnS35lXXF0KEtoS2llXXF1KEtoS39lXXF2KEtpS2plXXF3KEtpS4BlXXF4KEtqS2tlXXF5KEtqS4FlXXF6KEtrS4JlXXF7KEtMS01lXXF8KEtMS05lXXF9KEtyS1tlXXF+KEtfS2JlXXF/KEtjS2JlXXGAKEtkS2JlZVUFbGFiZWxxgUt8WAAAAAB9h1UIaGFsZmJvbmRxgkt8iH2HVQZyYWRpdXNxg0t8Rz/JmZmgAAAAfYdVC2xhYmVsT2Zmc2V0cYRLfE59h1UIZHJhd01vZGVxhUt8SwF9h1UIb3B0aW9uYWxxhn1VB2Rpc3BsYXlxh0t8SwJ9h3Uu'))
	crdInfo = cPickle.loads(base64.b64decode('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'))
	surfInfo = {'category': (0, None, {}), 'probeRadius': (0, None, {}), 'pointSize': (0, None, {}), 'name': [], 'density': (0, None, {}), 'colorMode': (0, None, {}), 'useLighting': (0, None, {}), 'transparencyBlendMode': (0, None, {}), 'molecule': [], 'smoothLines': (0, None, {}), 'lineWidth': (0, None, {}), 'allComponents': (0, None, {}), 'twoSidedLighting': (0, None, {}), 'customVisibility': [], 'drawMode': (0, None, {}), 'display': (0, None, {}), 'customColors': []}
	vrmlInfo = {'subid': (0, None, {}), 'display': (0, None, {}), 'id': (0, None, {}), 'vrmlString': [], 'name': (0, None, {})}
	colors = {u'': ((0.780392, 1, 0.780392), 1, u''), u'Ru': ((0.141176, 0.560784, 0.560784), 1, u'default'), u'Re': ((0.14902, 0.490196, 0.670588), 1, u'default'), u'Rf': ((0.8, 0, 0.34902), 1, u'default'), u'Ra': ((0, 0.490196, 0), 1, u'default'), u'Rb': ((0.439216, 0.180392, 0.690196), 1, u'default'), u'Rn': ((0.258824, 0.509804, 0.588235), 1, u'default'), u'Rh': ((0.0392157, 0.490196, 0.54902), 1, u'default'), u'Be': ((0.760784, 1, 0), 1, u'default'), u'Ba': ((0, 0.788235, 0), 1, u'default'), u'Bh': ((0.878431, 0, 0.219608), 1, u'default'), u'Bi': ((0.619608, 0.309804, 0.709804), 1, u'default'), u'Bk': ((0.541176, 0.309804, 0.890196), 1, u'default'), u'Br': ((0.65098, 0.160784, 0.160784), 1, u'default'), u'H': ((1, 1, 1), 1, u'default'), u'P': ((1, 0.501961, 0), 1, u'default'), u'Os': ((0.14902, 0.4, 0.588235), 1, u'default'), u'Ge': ((0.4, 0.560784, 0.560784), 1, u'default'), u'Gd': ((0.270588, 1, 0.780392), 1, u'default'), u'Ga': ((0.760784, 0.560784, 0.560784), 1, u'default'), u'Pr': ((0.85098, 1, 0.780392), 1, u'default'), u'Pt': ((0.815686, 0.815686, 0.878431), 1, u'default'),
u'Pu': ((0, 0.419608, 1), 1, u'default'), u'Mg': ((0.541176, 1, 0), 1, u'default'), u'Pb': ((0.341176, 0.34902, 0.380392), 1, u'default'), u'Pa': ((0, 0.631373, 1), 1, u'default'), u'Pd': ((0, 0.411765, 0.521569), 1, u'default'), u'Cd': ((1, 0.85098, 0.560784), 1, u'default'), u'Po': ((0.670588, 0.360784, 0), 1, u'default'), u'Pm': ((0.639216, 1, 0.780392), 1, u'default'), u'Hs': ((0.901961, 0, 0.180392), 1, u'default'), u'Ho': ((0, 1, 0.611765), 1, u'default'), u'Hf': ((0.301961, 0.760784, 1), 1, u'default'), u'Hg': ((0.721569, 0.721569, 0.815686), 1, u'default'), u'He': ((0.85098, 1, 1), 1, u'default'), u'Md': ((0.701961, 0.0509804, 0.65098), 1, u'default'), u'C': ((0.564706, 0.564706, 0.564706), 1, u'default'), u'K': ((0.560784, 0.25098, 0.831373), 1, u'default'), u'Mn': ((0.611765, 0.478431, 0.780392), 1, u'default'), u'O': ((1, 0.0509804, 0.0509804), 1, u'default'), u'Mt': ((0.921569, 0, 0.14902), 1, u'default'), u'S': ((1, 1, 0.188235), 1, u'default'), u'W': ((0.129412, 0.580392, 0.839216), 1, u'default'), u'sky blue': ((0.529412, 0.807843, 0.921569), 1, u'default'),
u'Zn': ((0.490196, 0.501961, 0.690196), 1, u'default'), u'plum': ((0.866667, 0.627451, 0.866667), 1, u'default'), u'Eu': ((0.380392, 1, 0.780392), 1, u'default'), u'Zr': ((0.580392, 0.878431, 0.878431), 1, u'default'), u'Er': ((0, 0.901961, 0.458824), 1, u'default'), u'Ni': ((0.313725, 0.815686, 0.313725), 1, u'default'), u'No': ((0.741176, 0.0509804, 0.529412), 1, u'default'), u'Na': ((0.670588, 0.360784, 0.94902), 1, u'default'), u'Nb': ((0.45098, 0.760784, 0.788235), 1, u'default'), u'Nd': ((0.780392, 1, 0.780392), 1, u'default'), u'Ne': ((0.701961, 0.890196, 0.960784), 1, u'default'), u'Np': ((0, 0.501961, 1), 1, u'default'), u'Fr': ((0.258824, 0, 0.4), 1, u'default'), u'Fe': ((0.878431, 0.4, 0.2), 1, u'default'), u'Fm': ((0.701961, 0.121569, 0.729412), 1, u'default'), u'B': ((1, 0.709804, 0.709804), 1, u'default'), u'F': ((0.564706, 0.878431, 0.313725), 1, u'default'), u'Sr': ((0, 1, 0), 1, u'default'), u'N': ((0.188235, 0.313725, 0.972549), 1, u'default'), u'Kr': ((0.360784, 0.721569, 0.819608), 1, u'default'), u'Si': ((0.941176, 0.784314, 0.627451), 1, u'default'),
u'Sn': ((0.4, 0.501961, 0.501961), 1, u'default'), u'Sm': ((0.560784, 1, 0.780392), 1, u'default'), u'V': ((0.65098, 0.65098, 0.670588), 1, u'default'), u'Sc': ((0.901961, 0.901961, 0.901961), 1, u'default'), u'Sb': ((0.619608, 0.388235, 0.709804), 1, u'default'), u'Sg': ((0.85098, 0, 0.270588), 1, u'default'), u'Se': ((1, 0.631373, 0), 1, u'default'), u'Co': ((0.941176, 0.564706, 0.627451), 1, u'default'), u'Cm': ((0.470588, 0.360784, 0.890196), 1, u'default'), u'Cl': ((0.121569, 0.941176, 0.121569), 1, u'default'), u'Ca': ((0.239216, 1, 0), 1, u'default'), u'Cf': ((0.631373, 0.211765, 0.831373), 1, u'default'), u'Ce': ((1, 1, 0.780392), 1, u'default'), u'Xe': ((0.258824, 0.619608, 0.690196), 1, u'default'), u'Lu': ((0, 0.670588, 0.141176), 1, u'default'), u'Cs': ((0.341176, 0.0901961, 0.560784), 1, u'default'), u'Cr': ((0.541176, 0.6, 0.780392), 1, u'default'), u'Cu': ((0.784314, 0.501961, 0.2), 1, u'default'), u'La': ((0.439216, 0.831373, 1), 1, u'default'), u'Li': ((0.8, 0.501961, 1), 1, u'default'), u'Tl': ((0.65098, 0.329412, 0.301961), 1, u'default'),
u'Tm': ((0, 0.831373, 0.321569), 1, u'default'), u'Lr': ((0.780392, 0, 0.4), 1, u'default'), u'Th': ((0, 0.729412, 1), 1, u'default'), u'Ti': ((0.74902, 0.760784, 0.780392), 1, u'default'), u'tan': ((0.823529, 0.705882, 0.54902), 1, u'default'), u'Te': ((0.831373, 0.478431, 0), 1, u'default'), u'Tb': ((0.188235, 1, 0.780392), 1, u'default'), u'Tc': ((0.231373, 0.619608, 0.619608), 1, u'default'), u'Ta': ((0.301961, 0.65098, 1), 1, u'default'), u'Yb': ((0, 0.74902, 0.219608), 1, u'default'), u'Db': ((0.819608, 0, 0.309804), 1, u'default'), u'Dy': ((0.121569, 1, 0.780392), 1, u'default'), u'At': ((0.458824, 0.309804, 0.270588), 1, u'default'), u'I': ((0.580392, 0, 0.580392), 1, u'default'), u'U': ((0, 0.560784, 1), 1, u'default'), u'Y': ((0.580392, 1, 1), 1, u'default'), u'Ac': ((0.439216, 0.670588, 0.980392), 1, u'default'), u'Ag': ((0.752941, 0.752941, 0.752941), 1, u'default'), u'Ir': ((0.0901961, 0.329412, 0.529412), 1, u'default'), u'Am': ((0.329412, 0.360784, 0.94902), 1, u'default'), u'Al': ((0.74902, 0.65098, 0.65098), 1, u'default'), u'As': ((0.741176, 0.501961, 0.890196), 1, u'default'),
u'Ar': ((0.501961, 0.819608, 0.890196), 1, u'default'), u'Au': ((1, 0.819608, 0.137255), 1, u'default'), u'Es': ((0.701961, 0.121569, 0.831373), 1, u'default'), u'In': ((0.65098, 0.458824, 0.45098), 1, u'default'), u'Mo': ((0.329412, 0.709804, 0.709804), 1, u'default')}
	materials = {u'default': ((0.85, 0.85, 0.85), 30), u'': ((0.85, 0.85, 0.85), 30)}
	pbInfo = {'category': [u'distance monitor'], 'bondInfo': [{'color': (0, None, {}), 'atoms': [], 'label': (0, None, {}), 'halfbond': (0, None, {}), 'labelColor': (0, None, {}), 'labelOffset': (0, None, {}), 'drawMode': (0, None, {}), 'display': (0, None, {})}], 'lineType': (1, 2, {}), 'color': (1, 58, {}), 'optional': {'fixedLabels': (True, False, (1, False, {}))}, 'display': (1, True, {}), 'showStubBonds': (1, False, {}), 'lineWidth': (1, 1, {}), 'stickScale': (1, 1, {}), 'id': [-2]}
	modelAssociations = {}
	colorInfo = (61, (u'', (0, 0, 0, 0.2)), {(u'', (0.193839, 0.452601, 0.218188, 1)): [5], (u'', (0.1184, 0.223266, 0.901463, 0.2)): [6, 9], (u'green', (0, 1, 0, 1)): [60], (u'', (0.663999, 0.164948, 0.651698, 1)): [4], (u'H', (1, 1, 1, 1)): [15], (u'N', (0.188235, 0.313725, 0.972549, 1)): [16], (u'', (1, 1, 1, 1)): [59], (u'S', (1, 1, 0.188235, 1)): [12], (u'O', (1, 0.0509804, 0.0509804, 1)): [13], (u'sky blue', (0.529412, 0.807843, 0.921569, 1)): [1], (u'', (0.193839, 0.452601, 0.218188, 0.2)): [8, 11], (u'tan', (0.823529, 0.705882, 0.54902, 1)): [0], (u'', (0.529549, 0.808011, 0.921828, 1)): [2], (u'Br', (0.65098, 0.160784, 0.160784, 1)): [14], (u'red', (1, 0, 0, 1)): [31], (u'yellow', (1, 1, 0, 1)): [58], (u'', (0.663999, 0.164948, 0.651698, 0.2)): [7, 10], (u'blue', (0, 0, 1, 1)): [30], (u'', (0.1184, 0.223266, 0.901463, 1)): [3]})
	viewerInfo = {'cameraAttrs': {'center': (46.179, 14.47249997139, 2.9510000367165), 'fieldOfView': 17.567064639105, 'nearFar': (15.820337462141, -12.771337686254), 'ortho': False, 'eyeSeparation': 50.8, 'focal': 2.9510000367165}, 'viewerAttrs': {'silhouetteColor': None, 'clipping': False, 'showSilhouette': False, 'showShadows': False, 'viewSize': 6.8744100148427, 'labelsOnTop': True, 'depthCueRange': (0.5, 1), 'silhouetteWidth': 2, 'singleLayerTransparency': True, 'shadowTextureSize': 2048, 'backgroundImage': [None, 1, 2, 1, 0, 0], 'backgroundGradient': [('Chimera default', [(1, 1, 1, 1), (0, 0, 1, 1)], 1), 1, 0, 0], 'depthCue': True, 'highlight': 0, 'scaleFactor': 0.75960428690837, 'angleDependentTransparency': True, 'backgroundMethod': 0}, 'viewerHL': 60, 'cameraMode': 'mono', 'detail': 1.5, 'viewerFog': None, 'viewerBG': 59}

	replyobj.status("Initializing session restore...", blankAfter=0,
		secondary=True)
	from SimpleSession.versions.v65 import expandSummary
	init(dict(enumerate(expandSummary(colorInfo))))
	replyobj.status("Restoring colors...", blankAfter=0,
		secondary=True)
	restoreColors(colors, materials)
	replyobj.status("Restoring molecules...", blankAfter=0,
		secondary=True)
	restoreMolecules(molInfo, resInfo, atomInfo, bondInfo, crdInfo)
	replyobj.status("Restoring surfaces...", blankAfter=0,
		secondary=True)
	restoreSurfaces(surfInfo)
	replyobj.status("Restoring VRML models...", blankAfter=0,
		secondary=True)
	restoreVRML(vrmlInfo)
	replyobj.status("Restoring pseudobond groups...", blankAfter=0,
		secondary=True)
	restorePseudoBondGroups(pbInfo)
	replyobj.status("Restoring model associations...", blankAfter=0,
		secondary=True)
	restoreModelAssociations(modelAssociations)
	replyobj.status("Restoring camera...", blankAfter=0,
		secondary=True)
	restoreViewer(viewerInfo)

try:
	restoreCoreModels()
except:
	reportRestoreError("Error restoring core models")

	replyobj.status("Restoring extension info...", blankAfter=0,
		secondary=True)


try:
	import StructMeasure
	from StructMeasure.DistMonitor import restoreDistances
	registerAfterModelsCB(restoreDistances, 1)
except:
	reportRestoreError("Error restoring distances in session")


def restoreMidasBase():
	formattedPositions = {'session-start': (1.0, 6.8744100148427, (46.179, 14.47249997139, 2.9510000367165), (15.820337462141, -12.771337686254), 2.9510000367165, {(0, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)), (102, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)), (2, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)), (101, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)), (1, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0)), (100, 0): ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))}, {(2, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0), (1, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0), (0, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0), (102, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0), (101, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0), (100, 0, 'Molecule'): (False, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False, 5.0)}, 4, (46.669999921321875, 9.811500097751617, 1.5244998879432687), False, 17.567064639105)}
	import Midas
	Midas.restoreMidasBase(formattedPositions)
try:
	restoreMidasBase()
except:
	reportRestoreError('Error restoring Midas base state')


def restoreMidasText():
	from Midas import midas_text
	midas_text.aliases = {}
	midas_text.userSurfCategories = {}

try:
	restoreMidasText()
except:
	reportRestoreError('Error restoring Midas text state')


def restore_volume_data():
 volume_data_state = \
  {
   'class': 'Volume_Manager_State',
   'data_and_regions_state': [ ],
   'version': 2,
  }
 from VolumeViewer import session
 session.restore_volume_data_state(volume_data_state)

try:
  restore_volume_data()
except:
  reportRestoreError('Error restoring volume data')


def restore_cap_attributes():
 cap_attributes = \
  {
   'cap_attributes': [ ],
   'cap_color': None,
   'cap_offset': 0.01,
   'class': 'Caps_State',
   'default_cap_offset': 0.01,
   'mesh_style': False,
   'shown': True,
   'subdivision_factor': 1.0,
   'version': 1,
  }
 import SurfaceCap.session
 SurfaceCap.session.restore_cap_attributes(cap_attributes)
registerAfterModelsCB(restore_cap_attributes)

geomData = {'AxisManager': {}, 'CentroidManager': {}, 'PlaneManager': {}}

try:
	from StructMeasure.Geometry import geomManager
	geomManager._restoreSession(geomData)
except:
	reportRestoreError("Error restoring geometry objects in session")


def restoreSession_RibbonStyleEditor():
	import SimpleSession
	import RibbonStyleEditor
	userScalings = []
	userXSections = []
	userResidueClasses = []
	residueData = [(6, 'Chimera default', 'rounded', u'unknown'), (7, 'Chimera default', 'rounded', u'unknown'), (8, 'Chimera default', 'rounded', u'unknown'), (9, 'Chimera default', 'rounded', u'unknown'), (10, 'Chimera default', 'rounded', u'unknown'), (11, 'Chimera default', 'rounded', u'unknown')]
	flags = RibbonStyleEditor.NucleicDefault1
	SimpleSession.registerAfterModelsCB(RibbonStyleEditor.restoreState,
				(userScalings, userXSections,
				userResidueClasses, residueData, flags))
try:
	restoreSession_RibbonStyleEditor()
except:
	reportRestoreError("Error restoring RibbonStyleEditor state")

trPickle = 'gAJjQW5pbWF0ZS5UcmFuc2l0aW9ucwpUcmFuc2l0aW9ucwpxASmBcQJ9cQMoVQxjdXN0b21fc2NlbmVxBGNBbmltYXRlLlRyYW5zaXRpb24KVHJhbnNpdGlvbgpxBSmBcQZ9cQcoVQZmcmFtZXNxCEsBVQ1kaXNjcmV0ZUZyYW1lcQlLAVUKcHJvcGVydGllc3EKXXELVQNhbGxxDGFVBG5hbWVxDVUMY3VzdG9tX3NjZW5lcQ5VBG1vZGVxD1UGbGluZWFycRB1YlUIa2V5ZnJhbWVxEWgFKYFxEn1xEyhoCEsUaAlLAWgKXXEUaAxhaA1VCGtleWZyYW1lcRVoD2gQdWJVBXNjZW5lcRZoBSmBcRd9cRgoaAhLAWgJSwFoCl1xGWgMYWgNVQVzY2VuZXEaaA9oEHVidWIu'
scPickle = 'gAJjQW5pbWF0ZS5TY2VuZXMKU2NlbmVzCnEBKYFxAn1xA1UHbWFwX2lkc3EEfXNiLg=='
kfPickle = 'gAJjQW5pbWF0ZS5LZXlmcmFtZXMKS2V5ZnJhbWVzCnEBKYFxAn1xA1UHZW50cmllc3EEXXEFc2Iu'
def restoreAnimation():
	'A method to unpickle and restore animation objects'
	# Scenes must be unpickled after restoring transitions, because each
	# scene links to a 'scene' transition. Likewise, keyframes must be 
	# unpickled after restoring scenes, because each keyframe links to a scene.
	# The unpickle process is left to the restore* functions, it's 
	# important that it doesn't happen prior to calling those functions.
	import SimpleSession
	from Animate.Session import restoreTransitions
	from Animate.Session import restoreScenes
	from Animate.Session import restoreKeyframes
	SimpleSession.registerAfterModelsCB(restoreTransitions, trPickle)
	SimpleSession.registerAfterModelsCB(restoreScenes, scPickle)
	SimpleSession.registerAfterModelsCB(restoreKeyframes, kfPickle)
try:
	restoreAnimation()
except:
	reportRestoreError('Error in Animate.Session')

def restoreLightController():
	import Lighting
	Lighting._setFromParams({'ratio': 1.25, 'brightness': 1.16, 'material': [30.0, (0.85, 0.85, 0.85), 1.0], 'back': [(0.3574067443365933, 0.6604015517481455, -0.6604015517481456), (1.0, 1.0, 1.0), 0.0], 'mode': 'two-point', 'key': [(-0.3574067443365933, 0.6604015517481455, 0.6604015517481456), (1.0, 1.0, 1.0), 1.0], 'contrast': 0.83, 'fill': [(0.2505628070857316, 0.2505628070857316, 0.9351131265310294), (1.0, 1.0, 1.0), 0.0]})
try:
	restoreLightController()
except:
	reportRestoreError("Error restoring lighting parameters")


def restoreRemainder():
	from SimpleSession.versions.v65 import restoreWindowSize, \
	     restoreOpenStates, restoreSelections, restoreFontInfo, \
	     restoreOpenModelsAttrs, restoreModelClip, restoreSilhouettes

	curSelIds =  []
	savedSels = []
	openModelsAttrs = { 'cofrMethod': 4 }
	windowSize = (904, 945)
	xformMap = {0: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True), 1: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True), 2: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True), 3: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True), 4: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True), 5: (((-0.18614921762158, 0.25844186963174, 0.94792207949816), 114.30009998269), (72.982139339236, -23.280282315404, 21.107007696256), True)}
	fontInfo = {'face': ('Fixed', 'Bold', 30)}
	clipPlaneInfo = {}
	silhouettes = {0: True, 1: True, 2: True, 3: True, 4: True, 5: True, 374: True}

	replyobj.status("Restoring window...", blankAfter=0,
		secondary=True)
	restoreWindowSize(windowSize)
	replyobj.status("Restoring open states...", blankAfter=0,
		secondary=True)
	restoreOpenStates(xformMap)
	replyobj.status("Restoring font info...", blankAfter=0,
		secondary=True)
	restoreFontInfo(fontInfo)
	replyobj.status("Restoring selections...", blankAfter=0,
		secondary=True)
	restoreSelections(curSelIds, savedSels)
	replyobj.status("Restoring openModel attributes...", blankAfter=0,
		secondary=True)
	restoreOpenModelsAttrs(openModelsAttrs)
	replyobj.status("Restoring model clipping...", blankAfter=0,
		secondary=True)
	restoreModelClip(clipPlaneInfo)
	replyobj.status("Restoring per-model silhouettes...", blankAfter=0,
		secondary=True)
	restoreSilhouettes(silhouettes)

	replyobj.status("Restoring remaining extension info...", blankAfter=0,
		secondary=True)
try:
	restoreRemainder()
except:
	reportRestoreError("Error restoring post-model state")
from SimpleSession.versions.v65 import makeAfterModelsCBs
makeAfterModelsCBs()

from SimpleSession.versions.v65 import endRestore
replyobj.status('Finishing restore...', blankAfter=0, secondary=True)
endRestore({})
replyobj.status('', secondary=True)
replyobj.status('Restore finished.')

