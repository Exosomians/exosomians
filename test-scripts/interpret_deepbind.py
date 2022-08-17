import exoNet
import numpy as np


data_name = "final"
motif_len = 10


sequences = np.load(f"./Data/{data_name}/sequences.npy", allow_pickle=True)[:, :, :4]
labels = np.load(f"./Data/{data_name}/labels.npy", allow_pickle=True).reshape((-1, ))

deepbind = exoNet.models.DeepBind(seq_len=data.shape[1],
                                  n_channels=data.shape[2],
                                  n_classes=2,
                                  use_batchnorm=True,
                                  lr=0.00001,
                                  padding="valid",
                                  model_path=f"./models/DeepBind/{data_name}/",
                                  dropout_rate=0.1,
                                 )

deepbind.train(data,
               labels, 
               le={"NO": 0, "YES": 1},
               n_epochs=1,
               batch_size=64,
               early_stopping_kwargs={"patience": 10, "monitor": "val_loss"},
              )


motif_scores = deepbind.aux_models['conv'].predict(data) # (None, seq_len - motif_len, n_filters=32)

position_frequency_matrix = np.zeros(shape=(32, motif_len, 4))
for i in range(32):
    for j in range(data.shape[0]):
        k = np.argmax(motif_scores[j, :, i]) # desired subsequence: [k:k+m]
        subsequence = sequences[j, [k:k+motif_len], :]
        position_frequency_matrix[i] += subsequence

np.save(arr=position_frequency_matrix, file="deepbind.pfm.npy")




