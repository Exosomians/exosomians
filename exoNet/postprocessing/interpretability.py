# from deepexplain.tensorflow import DeepExplain
import innvestigate
import innvestigate.utils as iutils
from keras import backend as K
from keras.utils import to_categorical
import exoNet

METHODS = {
    "DeepTaylor": {"name": "deep_taylor"},
    "DeepLIFT": {"name": "deeplift"}
}

def interpret(network, data, labels, method="DeepTaylor", **method_kwargs):
    method_dict = METHODS[method]
    method_name = method_dict.get("name", None)
    
    if method_name == "deeplift":
        with DeepExplain(session=K.get_session()) as de:
            input_tensor = network.model.layers[0].input
            fModel = network.model
            target_tensor = fModel(input_tensor)
            
            y_encoded, _ = exoNet.utils.label_encoder(labels, network.label_encoder)
            y_onehot = to_categorical(y_encoded, num_classes=2)
            
            attributions = de.explain('deeplift', target_tensor, input_tensor, data, ys=y_onehot, **method_kwargs)
            
            return attributions
    elif method_name == "deep_taylor":
        model_wo_sm = iutils.keras.graph.model_wo_softmax(network.model)
        analyzer = innvestigate.create_analyzer(method_name,        
                                                model_wo_sm, 
                                                )
        heatmaps = analyzer.analyze(data)
        return heatmaps
    else:
        raise Exception("Invalid interpretability method")