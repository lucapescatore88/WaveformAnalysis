import yaml

class objcfg :
    def __init__(self,entries) :
        self.__dict__.update(**entries)


def get_checked_cfg(fname) :

    print "Reading configuration from", fname
    with open(fname,"r") as cfgfile :
        cfg = yaml.load(cfgfile)
        for inc in cfg.get('includes',[]) :
            cfg.update(yaml.load(open(inc)))

        if len(cfg['VerticalOffset']) != len(cfg['Voltage']) :
            print 'VerticalOffset not same length as voltages'
            return
        if len(cfg['TriggerLevel']) != len(cfg['Voltage']) :
            print 'TriggerLevel not same length as voltages'
            return
        if len(cfg['VerticalScale']) != len(cfg['Voltage']) :
            print 'VerticalScale not same length as voltages'
            return
        print cfg 
     
    return objcfg(cfg)